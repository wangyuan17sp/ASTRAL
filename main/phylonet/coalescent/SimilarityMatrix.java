package phylonet.coalescent;


import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

/**
 * Implements a Distance method
 * @author smirarab
 *
 */
public class SimilarityMatrix extends AbstractMatrix implements Matrix {

    public SimilarityMatrix(int n) {
        this.n = n;
    }

    public SimilarityMatrix(float[][] from) {
        this.n = from.length;
        this.matrix = from;
    }

    int compareTwoValues(float f1, float f2) {
        int vc = Float.compare(f1,f2);
        return - vc;
    }
    

    public int getBetterSideByFourPoint(int x, int a, int b, int c) {
        double xa = this.matrix[x][a];
        double xb = this.matrix[x][b];
        double xc = this.matrix[x][c];
        double ab = this.matrix[a][b];
        double ac = this.matrix[a][c];
        double bc = this.matrix[b][c];
        double ascore = xa + bc  - (xb + ac); // Note this is similartiy, not distance
        double bscore = xb + ac  - (xa + bc); 
        double cscore = xc + ab - (xb + ac); 
        return ascore >= bscore ?
                ascore >= cscore ? a : c :
                    bscore >= cscore ? b : c;	
    }


    List<BitSet> UPGMA() {
        
        List<BitSet> bsList = new ArrayList<BitSet>(n);
        List<TreeSet<Integer>> indsBySim = new ArrayList<TreeSet<Integer>>(n);
        List<float[]> sims = new ArrayList<float[]>(n);
        List<Integer> range = Utils.getRange(n);
        List<Integer> weights = Utils.getOnes(n);
    
        for (int i = 0; i< n; i++) {
            BitSet bs = new BitSet();
            bs.set(i);
            bsList.add(bs);
            final float[] is = this.matrix[i].clone();
    
            sims.add(is);
            range.remove(i);
            TreeSet<Integer> sortColumn = this.sortColumn(range, is);
            range.add(i,i);
            indsBySim.add(sortColumn);
        }
        
    
        return upgmaLoop(weights, bsList, indsBySim, sims, n, false);
    }

    List<BitSet> resolveByUPGMA(List<BitSet> bsList, boolean original) {
    
        List<BitSet> internalBSList;
        if (original) {
            internalBSList = new ArrayList<BitSet>(bsList);
        } else {
            internalBSList = new ArrayList<BitSet>();
        }
    
        int size = bsList .size();
        List<TreeSet<Integer>> indsBySim = new ArrayList<TreeSet<Integer>>(size);
        List<float[]> sims = new ArrayList<float[]>(size);
        List<Integer> range = Utils.getRange(size);
        List<Integer> weights = new ArrayList<Integer>(size);
    
        for (int i = 0; i < size; i++) {
            if (!original) {
                BitSet internalBS = new BitSet(size);
                internalBS.set(i);
                internalBSList.add(internalBS);
            }
    
            final float[] is = new float[size];// this.similarityMatrix[i].clone();
            BitSet bsI = bsList.get(i);
    
            weights.add(bsI.cardinality());
            sims.add(is);
    
            for (int j = 0; j < size; j++) {
    
                BitSet bsJ = bsList.get(j);
                int c = 0;
                if (i == j) {
                    is[j] = 1;
                    continue;
                }
                for (int k = bsI.nextSetBit(0); k >= 0; k = bsI.nextSetBit(k + 1)) {
                    for (int l = bsJ.nextSetBit(0); l >= 0; l = bsJ.nextSetBit(l + 1)) {
                        //						System.err.println("k :"+k+" l : "+l);
                        is[j] += this.matrix[k][l];
                        c++;
                    }
                }
                if (c == 0) {
                    throw new RuntimeException("Error: "+bsI + " "+bsJ);
                }
                is[j] /= c;
            }
    
            range.remove(i);
            TreeSet<Integer> sortColumn = this.sortColumn(range, is);
            range.add(i,i);
            indsBySim.add(sortColumn);
        }
    
        return upgmaLoop(weights, internalBSList, indsBySim, sims, size, false);
    }

    void populateByQuartetDistance(List<STITreeCluster> treeAllClusters, List<Tree> geneTrees) {
    
        this.matrix = new float[n][n];
        long [][] denom = new long [n][n];
    
        int k = 0;
        for (Tree tree :  geneTrees) {
    
            for (TNode node : tree.postTraverse()) {
                if (node.isLeaf()) {
                    BitSet tmp = new BitSet(n);
                    tmp.set(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
                    ((STINode)node).setData(tmp);
                } else {
    
                    BitSet newbs = new BitSet(n);
                    for (TNode cn: node.getChildren()) {
                        BitSet c = (BitSet) ((STINode)cn).getData();
                        newbs.or(c);
                    }
    
                    ((STINode)node).setData(newbs);
    
                }
            }
        }
    
        for (Tree tree :  geneTrees) {
            STITreeCluster treeallCL = treeAllClusters.get(k++);
    
            Integer treeall = treeallCL.getClusterSize();
    
            for (TNode node : tree.postTraverse()) {
                if (node.isLeaf()) { 
                    continue;
                }
                BitSet cluster = (BitSet) ((STINode)node).getData();
                BitSet others = (BitSet) treeallCL.getBitSet().clone();
                others.andNot(cluster);
                ArrayList<BitSet> children = new ArrayList<BitSet>();
                long totalPairs = 0;
                long totalUnresolvedPairs = 0;
                for (TNode cn: node.getChildren()) {
                    BitSet c = (BitSet) ((STINode)cn).getData();
                    children.add(c);
                    long cc = c.cardinality();
                    totalPairs += cc*(cc-1);
                    totalUnresolvedPairs += cc * (treeall - cc); 
                }
                if (others.cardinality() != 0) {
                    children.add(others);
                    long cc = others.cardinality();
                    totalPairs += cc*(cc-1);
                    totalUnresolvedPairs += cc * (treeall - cc);
                }
                totalPairs /= 2;
                totalUnresolvedPairs /= 2;
    
    
                for (int j = 0; j < children.size(); j++ ) {
                    BitSet left = children.get(j);
                    long lc = left.cardinality();
                    long lcu = lc * (treeall - lc);
                    long lcp = lc*(lc-1)/2;
                    for (int i = j+1; i < children.size(); i++ ) {
                        BitSet right = children.get(i);
                        long rc = right.cardinality();
                        long rcu = rc * (treeall - lc - rc);
                        long rcp = rc*(rc-1)/2;
                        double sim = (totalPairs - lcp - rcp) // the number of fully resolved quartets
                                //+ (totalUnresolvedPairs - lcu - rcu) / 3.0 // we count partially resolved quartets
                                ; 
                        updateQuartetDistanceTri( left, right, matrix, sim);
                    }
                }
            }
    
    
            BitSet all = treeallCL.getBitSet();
            int c = all.cardinality() - 2;
            for (int l = all.nextSetBit(0); l >= 0; l=all.nextSetBit(l+1)) {
                for (int r = all.nextSetBit(0); r >= 0; r=all.nextSetBit(r+1)) {
                    denom[l][r] += c*(c-1)/2;
                    denom[r][l] = denom[l][r];
                }
            }
    
        }
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                if (denom[i][j] == 0)
                    matrix[i][j] = 0;
                else
                    matrix[i][j] = matrix[i][j] / (denom[i][j]/2);
                if (i == j) {
                    matrix[i][j] = 1;
                }
                matrix[j][i] = matrix[i][j];
            }
            //System.err.println(Arrays.toString(similarityMatrix[i]));
        }
    }

    @Override
    public boolean isDistance() {
        return false;
    }

    @Override
    public List<BitSet> inferTreeBitsets() {
        
        return UPGMA();
    }

    @Override
    public Matrix populate(List<STITreeCluster> treeAllClusters, List<Tree> geneTrees, SpeciesMapper spm) {
        this.populateByQuartetDistance(treeAllClusters, geneTrees);
        return spm.convertToSpeciesDistance(this);
    }

    @Override
    public List<BitSet> resolvePolytomy(List<BitSet> bsList, boolean original) {
        return resolveByUPGMA(bsList, original);
    }

    @Override
    Matrix factory(float[][] from) {
        return new SimilarityMatrix(from);
    }

    private List<BitSet> upgmaLoop(List<Integer> weights, List<BitSet> bsList,
            List<TreeSet<Integer>> indsBySim, List<float[]> sims, int left, boolean randomize) {
    
        List<BitSet> ret = new ArrayList<BitSet>();
        while ( left > 2) {
            int closestI = -1;
            int closestJ = -1;
            float bestHit = -1;
            for (int i = 0; i < indsBySim.size(); i++) {
                if (indsBySim.get(i) == null)
                    continue;
                int j = indsBySim.get(i).first();
                if (sims.get(i)[j] > bestHit || (randomize & sims.get(i)[i] == bestHit & GlobalMaps.random.nextBoolean())) {
                    bestHit = sims.get(i)[j];
                    closestI = i;
                    closestJ = j;
                }
            }
            BitSet bs = (BitSet) bsList.get(closestI).clone();
            bs.or(bsList.get(closestJ));
            bsList.set(closestJ,null);
            bsList.set(closestI,bs);
    
            float[] jDist = sims.get(closestJ);
            float[] iDist = sims.get(closestI).clone();
            for (int k = 0; k < sims.size(); k++) {
                if (k == closestJ || sims.get(k) == null) {
                    continue;
                }
    
                if ( k != closestI) {
                    float newSimToI = (iDist[k] * weights.get(closestI) + jDist[k] * weights.get(closestJ)) /
                            ( weights.get(closestI) + weights.get(closestJ));
    
                    indsBySim.get(k).remove(closestI);
                    sims.get(k)[closestI] = newSimToI;
                    indsBySim.get(k).add(closestI);
    
                    indsBySim.get(closestI).remove(k);
                    sims.get(closestI)[k] = newSimToI;
                    indsBySim.get(closestI).add(k);
                }
    
                indsBySim.get(k).remove(closestJ);
                sims.get(k)[closestJ] = -1;
                //indsBySim.get(k).add(closestJ);
            }
    
            sims.set(closestJ,null);
            indsBySim.set(closestJ,null);
            weights.set(closestI, weights.get(closestI) + weights.get(closestJ));
            weights.set(closestJ,null);
            ret.add(bs);
            left--;
        }
        return ret;
    }

    /*
    private void updateQuartetDistanceForPair (Integer treeall, BitSet left,
    		BitSet right, float[][] matrix) {
    	long c = treeall - left.cardinality() - right.cardinality();
    	c = c*(c-1)/2;
    	for (int l = left.nextSetBit(0); l >= 0; l=left.nextSetBit(l+1)) {
    		for (int r = right.nextSetBit(0); r >= 0; r=right.nextSetBit(r+1)) {
    			matrix[l][r] += c;
    			matrix[r][l] = matrix[l][r];
    		}
    	}
    }
     */
    private void updateQuartetDistanceTri(BitSet left,
            BitSet right, float[][] matrix,double d) {
        if (d == 0)
            return;
        for (int l = left.nextSetBit(0); l >= 0; l=left.nextSetBit(l+1)) {
            for (int r = right.nextSetBit(0); r >= 0; r=right.nextSetBit(r+1)) {
                matrix[l][r] += d;
                matrix[r][l] = matrix[l][r];
            }
        }
    }

    
    /*
    public static void main(String[] args) {
        String newick;
        try {
            File test_tree = new File("/Users/Yuan/Desktop/UCSD-Summer17/ASTRAL/ASTRAL/main/test_data/test.tt");
            BufferedReader in = new BufferedReader(new FileReader(test_tree));
            newick = in.readLine();
            System.out.println(newick);
            in.close();
        } catch (IOException e) { throw new RuntimeException("Cannot find file: " + e); }

        Tree true_tree = null;

        try {
            NewickReader nr = new NewickReader(new StringReader(newick));
            true_tree = nr.readTree();
            System.out.println(true_tree);
            String[] leaves = true_tree.getLeaves();
            for (int i = 0; i < leaves.length; i++) {
                GlobalMaps.taxonIdentifier.taxonId(leaves[i]);
            }
        } catch (ParseException e) {
            throw new RuntimeException("Failed to Parse Tree: " , e);
        } catch (IOException e) {
            throw new RuntimeException();
        }

        List<Tree> tree = new ArrayList<Tree>();
        tree.add(true_tree);
        System.out.println(tree);
        SimilarityMatrix geneMatrix = new SimilarityMatrix(GlobalMaps.taxonIdentifier.taxonCount());
        System.out.println(GlobalMaps.taxonIdentifier.taxonCount());
        GlobalMaps.taxonNameMap = new TaxonNameMap();
        SpeciesMapper spm = GlobalMaps.taxonNameMap.getSpeciesIdMapper();
        SimilarityMatrix speciesMatrix = geneMatrix.matricesByBranchDistance(tree, spm);

        String newick_p;
        try {
            File test_tree_p = new File("/Users/Yuan/Desktop/UCSD-Summer17/ASTRAL/ASTRAL/main/test_data/poly.tt");
            BufferedReader in = new BufferedReader(new FileReader(test_tree_p));
            newick = in.readLine();
            in.close();
        } catch (IOException e) { throw new RuntimeException("Cannot find file: " + e); }

        Tree poly_tree = null;

        try {
            NewickReader nr_p = new NewickReader(new StringReader(newick));
            poly_tree = nr_p.readTree();

        } catch (ParseException e) {
            throw new RuntimeException("Failed to Parse Tree: " , e);
        } catch (IOException e) {
            throw new RuntimeException();
        }

        List<BitSet> bs1 = new ArrayList<BitSet>();
        List<BitSet> bs3 = new ArrayList<BitSet>();

        for (TNode node : poly_tree.postTraverse()) {
            if (node.isRoot()) {
                continue;
            } else if (node.isLeaf()) {
                BitSet leaf = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
                leaf.set(GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTTaxonIdentifier().taxonId(node.getName()));
                ((STINode)node).setData(leaf);
                System.out.println("Leaf is: " + leaf);
//                bs1.add(leaf);
            } else {
                System.out.println("OK! I am an internal node.");
                BitSet newbs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());

                for (TNode cn : node.getChildren()) {
                    BitSet c = (BitSet) ((STINode)cn).getData();
                    newbs.or(c);
                    if (node.getChildCount() > 2) {
                        bs3.add(c);
                    }
                }
                ((STINode)node).setData(newbs);
                bs1.add(newbs);
            }
        }
        System.out.println("BitSet1 is : " + bs1);
        System.out.println("Polytomy tree is: " + poly_tree);
        for (int i = 0; i < speciesMatrix.speciesDM.length; i++) {
            for (int j = 0; j < speciesMatrix.speciesDM.length; j++) {
                System.out.print(speciesMatrix.speciesDM[i][j] + " ");
            }
            System.out.println();
        }

        List<BitSet> bs2 = speciesMatrix.resolveByPhyDstar(bs1, true);
        bs2.addAll(bs3);
        System.out.println("BitSet2 is: " + bs2);
        List<STITreeCluster> phyd = new ArrayList<STITreeCluster>();
        for(int i = 0; i < bs2.size(); i++) {
            STITreeCluster sti = new STITreeCluster(GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTTaxonIdentifier());
            sti.setCluster(bs2.get(i));
            phyd.add(sti);
        }
        Tree result = Utils.buildTreeFromClusters(phyd, GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTTaxonIdentifier(), false);
        System.out.println("I am the final result: " + result.toNewick());

    } */
}
