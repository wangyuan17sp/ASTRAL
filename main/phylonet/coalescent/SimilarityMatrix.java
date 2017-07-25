package phylonet.coalescent;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeSet;

import phylonet.tree.io.NewickReader;
import phylonet.tree.io.ParseException;
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
public class SimilarityMatrix {

    private float[][] similarityMatrix;
    private float[][] speciesDM;
    private List<TreeSet<Integer>> orderedTaxonBySimilarity;
    private Integer n;

    public SimilarityMatrix(int n) {
        this.n = n;
    }

    public SimilarityMatrix(float[][] from) {
        this.n = from.length;
        this.similarityMatrix = from;
    }

    public int getSize() {
        return n;
    }

    public float get(int i, int j) {
        return this.similarityMatrix[i][j];
    }

    int getBetterSideByFourPoint(int x, int a, int b, int c) {
        double xa = this.similarityMatrix[x][a];
        double xb = this.similarityMatrix[x][b];
        double xc = this.similarityMatrix[x][c];
        double ab = this.similarityMatrix[a][b];
        double ac = this.similarityMatrix[a][c];
        double bc = this.similarityMatrix[b][c];
        double ascore = xa + bc  - (xb + ac); // Note this is similartiy, not distance
        double bscore = xb + ac  - (xa + bc); 
        double cscore = xc + ab - (xb + ac); 
        return ascore >= bscore ?
                ascore >= cscore ? a : c :
                    bscore >= cscore ? b : c;	
    }

    private List<TreeSet<Integer>> sortByDistance(float[][] refMatrix) {
        List<TreeSet<Integer>> ret = new ArrayList<TreeSet<Integer>>(n);
        List<Integer> range = Utils.getRange(n);
        for (int i = 0; i < n; i++) {
            final float[] js = refMatrix[i];
            TreeSet<Integer> indices = sortColumn(range, js);
            ret.add(indices);
        }
        return ret;
    }

    private TreeSet<Integer> sortColumn(List<Integer> range, final float[] js) {
        TreeSet<Integer> indices = new TreeSet<Integer>(new Comparator<Integer>() {

            @Override
            public int compare(Integer o1, Integer o2) {
                if (o1 == o2) {
                    return 0;
                }
                int comp = Float.compare(js[o1], js[o2]) ;
                return  comp == 0 ? - o1.compareTo(o2) : - comp;
            }
        });
        indices.addAll(range);
        return indices;
    }

    private void assureOrderedTaxa () {
        if (this.orderedTaxonBySimilarity == null) {
            this.orderedTaxonBySimilarity = this.sortByDistance(this.similarityMatrix);
        }
    }

    /**
     * Returns the id of the closest taxon that is either present in presentBS or
     * has a smaller id than mssingId
     * @param presentBS
     * @param missingId
     * @return
     */
    int getClosestPresentTaxonId(BitSet presentBS, int missingId) {
        this.assureOrderedTaxa();
        int closestId = -1;
        for (Integer other: this.orderedTaxonBySimilarity.get(missingId)){
            if ( missingId > other // other is already added
                    || presentBS.get(other) // other was in original tree
                    ) {
                closestId = other;
                break;
            }
        }
        if (closestId == -1) {
            throw new RuntimeException("Bug: this should not be reached");
        }
        return closestId;
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

    // 339023690, 11085
    void pupulateByBranchDistanceBS(List<Tree> geneTrees) {
        this.similarityMatrix = new float[n][n];
        float[][] pairNumMatrix = new float[n][n];

        for (Tree tree : geneTrees) {
            boolean[][] indicatorMatrix = new boolean[n][n];
            HashMap<Integer, Integer> distanceMap = new HashMap<Integer, Integer>();

            for (TNode node : tree.postTraverse()) {
                if (node.isLeaf()) {

                    // Setup BitSet
                    BitSet tmp = new BitSet(n);
                    tmp.set(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
                    ((STINode)node).setData(tmp);

                    // Update Map
                    distanceMap.put(GlobalMaps.taxonIdentifier.taxonId(node.getName()), 0);
                } else {

                    // Setup Bitset
                    List<BitSet> children = new ArrayList<BitSet>();
                    BitSet newbs = new BitSet(n);
                    for (TNode cn: node.getChildren()) {
                        BitSet c = (BitSet) ((STINode)cn).getData();
                        children.add(c);
                        newbs.or(c);
                    }
                    ((STINode)node).setData(newbs);

                    // Update similarity matrix
                    for (int i = 0; i < children.size() - 1; i++) {
                        BitSet left = children.get(i);
                        for (int j = i + 1; j < children.size(); j++) {
                            BitSet right = children.get(j);
                            for (int k = 0; k < left.length(); k++) {
                                if (left.get(k)) {
                                    for (int l = 0; l < right.length(); l++) {
                                        if (right.get(l)) {
                                            similarityMatrix[k][l] += distanceMap.get(k) + distanceMap.get(l) + 2;
                                            similarityMatrix[l][k] = similarityMatrix[k][l];
                                            pairNumMatrix[k][l] += 1;
                                            pairNumMatrix[l][k] = pairNumMatrix[k][l];
                                        }
                                    }
                                }
                            }
                        }
                    }


                    for (int index = 0; index < newbs.length(); index++) {
                        if (newbs.get(index)) {
                            distanceMap.put(index, distanceMap.get(index) + 1);
                        }
                    }
                }
            }	        
        }

        System.out.println(GlobalMaps.taxonIdentifier.taxonCount());
        for (int i = 0; i < n; i++) {
            System.out.print(GlobalMaps.taxonIdentifier.getTaxonName(i) + " ");
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    similarityMatrix[i][j] = 0;
                } else {
                    if (pairNumMatrix[i][j] != 0) {
                        similarityMatrix[i][j] /= pairNumMatrix[i][j];
                        similarityMatrix[i][j] = similarityMatrix[i][j];
                    } else {
                        similarityMatrix[i][j] = -99; 
                    }
                }
                System.out.print(similarityMatrix[i][j] + " ");         
            }
            System.out.println();
        }
    }

    void pupulateByBranchDistanceBL(List<Tree> geneTrees) {
        this.similarityMatrix = new float[n][n];
        float[][] pairNumMatrix = new float[n][n];

        for (Tree tree : geneTrees) {
            boolean[][] indicatorMatrix = new boolean[n][n];
            HashMap<Integer, Integer> distanceMap = new HashMap<Integer, Integer>();

            for (TNode node : tree.postTraverse()) {
                if (node.isLeaf()) {
                    distanceMap.put(GlobalMaps.taxonIdentifier.taxonId(node.getName()), 0);
                } else {
                    for (TNode cn1 : node.getLeaves()) {
                        for (TNode cn2 : node.getLeaves()) {
                            int cn1ID = GlobalMaps.taxonIdentifier.taxonId(cn1.getName());
                            int cn2ID = GlobalMaps.taxonIdentifier.taxonId(cn2.getName());
                            if (cn1ID == cn2ID || indicatorMatrix[cn1ID][cn2ID]) {
                                continue;
                            }
                            similarityMatrix[cn1ID][cn2ID] += distanceMap.get(cn1ID) + distanceMap.get(cn2ID) + 2;
                            indicatorMatrix[cn1ID][cn2ID] = true;
                            pairNumMatrix[cn1ID][cn2ID] += 1;
                        }
                    }

                    for (TNode cn : node.getLeaves()) {
                        int cnID = GlobalMaps.taxonIdentifier.taxonId(cn.getName());
                        distanceMap.put(cnID, distanceMap.get(cnID) + 1);
                    }

                }
            }           
        }


        for (int i = 0; i < n; i++) {
            System.err.print(GlobalMaps.taxonIdentifier.getTaxonName(i) + " ");
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    similarityMatrix[i][j] = 0;
                } else {
                    if (pairNumMatrix[i][j] != 0) {
                        similarityMatrix[i][j] /= pairNumMatrix[i][j];
                        similarityMatrix[i][j] = n - similarityMatrix[i][j];
                    } else {
                        similarityMatrix[i][j] = -1; 
                    }
                }
                System.err.print(similarityMatrix[i][j] + " ");         
            }
            System.err.println();
        }
    }


    void populateByQuartetDistance(List<STITreeCluster> treeAllClusters, List<Tree> geneTrees) {

        this.similarityMatrix = new float[n][n];
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
                        updateQuartetDistanceTri( left, right, similarityMatrix, sim);
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
                    similarityMatrix[i][j] = 0;
                else
                    similarityMatrix[i][j] = similarityMatrix[i][j] / (denom[i][j]/2);
                if (i == j) {
                    similarityMatrix[i][j] = 1;
                }
                similarityMatrix[j][i] = similarityMatrix[i][j];
            }
            //System.err.println(Arrays.toString(similarityMatrix[i]));
        }
    }

    // Not consider missing data, 11321
    SimilarityMatrix matricesByBranchDistance(List<Tree> geneTrees, SpeciesMapper spm) {
        this.similarityMatrix = new float[n][n]; // master gene matrix
        float[][] speciesSimilarityMatrix = new float[spm.getSpeciesCount()][spm.getSpeciesCount()]; // master species matrix
        float[][] speciesDistanceMatrix = new float[spm.getSpeciesCount()][spm.getSpeciesCount()];
        float[][] pairNumMatrix = new float[n][n];

        for (Tree tree : geneTrees) {
            HashMap<Integer, Integer> distanceMap = new HashMap<Integer, Integer>();
            float[][] geneSimilarityMatrix = new float[n][n];
            float[][] geneDissimilarityM = new float[n][n];
            float max = 0;

            for (TNode node : tree.postTraverse()) {
                if (node.isLeaf()) {
                    BitSet tmp = new BitSet(n);
                    tmp.set(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
                    ((STINode)node).setData(tmp);

                    distanceMap.put(GlobalMaps.taxonIdentifier.taxonId(node.getName()), 0);
                } else {
                    List<BitSet> children = new ArrayList<BitSet>();
                    BitSet newbs = new BitSet(n);
                    for (TNode cn: node.getChildren()) {
                        BitSet c = (BitSet) ((STINode)cn).getData();
                        children.add(c);
                        newbs.or(c);
                    }
                    ((STINode)node).setData(newbs);

                    for (int i = 0; i < children.size() - 1; i++) {
                        BitSet left = children.get(i);
                        for (int j = i + 1; j < children.size(); j++) {
                            BitSet right = children.get(j);
                            for (int k = 0; k < left.length(); k++) {
                                if (left.get(k)) {
                                    for (int l = 0; l < right.length(); l++) {
                                        if (right.get(l)) {
                                            geneDissimilarityM[k][l] = distanceMap.get(k) + distanceMap.get(l) + 2;
                                            geneDissimilarityM[l][k] = geneDissimilarityM[k][l];
                                            this.similarityMatrix[k][l] += geneDissimilarityM[k][l];
                                            this.similarityMatrix[l][k] = this.similarityMatrix[k][l];
                                            pairNumMatrix[k][l] += 1;
                                            pairNumMatrix[l][k] = pairNumMatrix[k][l];
                                            geneSimilarityMatrix[k][l] = n - geneDissimilarityM[k][l];
                                            geneSimilarityMatrix[l][k] = geneSimilarityMatrix[k][l];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    for (int index = 0; index < newbs.length(); index++) {
                        if (newbs.get(index)) {
                            distanceMap.put(index, distanceMap.get(index) + 1);
                        }
                    }
                }
            }
           
            SimilarityMatrix tmp = convertToSpeciesDistance2(spm, geneSimilarityMatrix);
            SimilarityMatrix speciesDM = convertToSpeciesDistance2(spm, geneDissimilarityM);
            for (int i = 0; i < spm.getSpeciesCount(); i++) {
                for (int j = i; j < spm.getSpeciesCount(); j++) {
                    speciesSimilarityMatrix[i][j] += tmp.similarityMatrix[i][j];
                    speciesSimilarityMatrix[j][i] = speciesSimilarityMatrix[i][j];
                    speciesDistanceMatrix[i][j] += speciesDM.similarityMatrix[i][j];
                    speciesDistanceMatrix[j][i] = speciesDistanceMatrix[i][j];
                }
            }
        }

//        System.out.println(GlobalMaps.taxonIdentifier.taxonCount());
        for (int i = 0; i < n; i++) {
//            System.out.print(GlobalMaps.taxonIdentifier.getTaxonName(i) + " ");
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    this.similarityMatrix[i][j] = 0;
                } else {
                    if (pairNumMatrix[i][j] != 0) {
                        this.similarityMatrix[i][j] /= pairNumMatrix[i][j];
                        this.similarityMatrix[i][j] = n - this.similarityMatrix[i][j];
                    } else {
                        this.similarityMatrix[i][j] = -99;
                    }
                } 
//                System.out.print(similarityMatrix[i][j] + " ");
            }
//            System.out.println();
        }

//        System.out.println(spm.getSpeciesCount());
        for (int i = 0; i < spm.getSpeciesCount(); i++) {
//            System.out.print(spm.getSpeciesName(i) + " ");
            for (int j = 0; j < spm.getSpeciesCount(); j++) {
                if (i == j) {
                    speciesSimilarityMatrix[i][j] = 0;
                    speciesDistanceMatrix[i][j] = 0;
                } else {
                    if (pairNumMatrix[i][j] != 0) {
                        speciesSimilarityMatrix[i][j] /= pairNumMatrix[i][j];
                        speciesDistanceMatrix[i][j] /= pairNumMatrix[i][j];
                    } else {
                        speciesSimilarityMatrix[i][j] = -99;
                        speciesDistanceMatrix[i][j] = -99;
                    }
                } 
//                System.out.print(speciesDistanceMatrix[i][j] + " ");
            }
//            System.out.println();
        }

        SimilarityMatrix ret = new SimilarityMatrix(speciesSimilarityMatrix.length);
        ret.similarityMatrix = speciesSimilarityMatrix;
        ret.speciesDM = speciesDistanceMatrix;
        this.speciesDM = speciesDistanceMatrix;
        
        return ret;
    }

    SimilarityMatrix convertToSpeciesDistance2(SpeciesMapper spm, float[][] similarityMatrix) {
        float [][] STsimMatrix = new float[spm.getSpeciesCount()][spm.getSpeciesCount()];
        float[][] denum = new float[spm.getSpeciesCount()][spm.getSpeciesCount()];
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                int stI =  spm.getSpeciesIdForTaxon(i);
                int stJ =  spm.getSpeciesIdForTaxon(j);
                STsimMatrix[stI][stJ] += similarityMatrix[i][j]; 
                STsimMatrix[stJ][stI] = STsimMatrix[stI][stJ];
                denum[stI][stJ] ++;
                denum[stJ][stI] ++;
            }
        }

        for (int i = 0; i < spm.getSpeciesCount(); i++) {
            for (int j = 0; j < spm.getSpeciesCount(); j++) {
                STsimMatrix[i][j] = denum[i][j] == 0 ? 0 : 
                    STsimMatrix[i][j] / denum[i][j];
            }
            STsimMatrix[i][i] = 0;
            //System.err.println(Arrays.toString(this.distSTMatrix[i]));
        }
        SimilarityMatrix ret = new SimilarityMatrix(STsimMatrix.length);
        ret.similarityMatrix = STsimMatrix;

        return ret;
    }


    SimilarityMatrix convertToSpeciesDistance(SpeciesMapper spm) {
        float [][] STsimMatrix = new float[spm.getSpeciesCount()][spm.getSpeciesCount()];
        float[][] denum = new float[spm.getSpeciesCount()][spm.getSpeciesCount()];
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                int stI =  spm.getSpeciesIdForTaxon(i);
                int stJ =  spm.getSpeciesIdForTaxon(j);
                STsimMatrix[stI][stJ] += this.similarityMatrix[i][j]; 
                STsimMatrix[stJ][stI] = STsimMatrix[stI][stJ];
                denum[stI][stJ] ++;
                denum[stJ][stI] ++;
            }
        }

        for (int i = 0; i < spm.getSpeciesCount(); i++) {
            for (int j = 0; j < spm.getSpeciesCount(); j++) {
                STsimMatrix[i][j] = denum[i][j] == 0 ? -1 : 
                    STsimMatrix[i][j] / denum[i][j];
            }
            STsimMatrix[i][i] = 0;
            //System.err.println(Arrays.toString(this.distSTMatrix[i]));
        }
        System.err.println("Species tree distances calculated ...");

        SimilarityMatrix ret = new SimilarityMatrix(STsimMatrix.length);
        ret.similarityMatrix = STsimMatrix;

        return ret;
    }
    SimilarityMatrix getInducedMatrix(HashMap<String, Integer> randomSample, TaxonIdentifier id) {

        int sampleSize = randomSample.size();
        float[][] sampleSimMatrix = new float[sampleSize][sampleSize];

        for (Entry<String, Integer> row : randomSample.entrySet()) {
            int rowI = id.taxonId(row.getKey());
            int i = row.getValue();
            for (Entry<String, Integer> col : randomSample.entrySet()) {
                int colJ = id.taxonId(col.getKey());
                sampleSimMatrix[i][col.getValue()] = this.similarityMatrix[rowI][colJ];
            }
        }
        SimilarityMatrix ret = new SimilarityMatrix(sampleSize);
        ret.similarityMatrix = sampleSimMatrix;
        return ret;
    }

    /*
     * 
     * TODO check what is this
     */

    SimilarityMatrix getInducedMatrix(List<Integer> sampleOrigIDs) {

        int sampleSize = sampleOrigIDs.size();
        SimilarityMatrix ret = new SimilarityMatrix(sampleSize);
        ret.similarityMatrix = new float [sampleSize][sampleSize];

        int i = 0;
        for (Integer rowI : sampleOrigIDs) {
            int j = 0;
            for (Integer colJ : sampleOrigIDs) {
                ret.similarityMatrix[i][j] = this.similarityMatrix[rowI][colJ];
                j++;
            }
            i++;
        }
        return ret;
    }

    //TODO: generate iterable, not list
    Iterable<BitSet> getQuadraticBitsets() {
        List<BitSet> newBitSets = new ArrayList<BitSet>();
        ArrayList<Integer> inds = new ArrayList<Integer> (n);
        for (int i = 0; i < n; i++) {
            inds.add(i);
        }
        for (final float[] fs : this.similarityMatrix) {
            Collections.sort(inds, new Comparator<Integer>() {

                @Override
                public int compare(Integer i1, Integer i2) {
                    if (i1 == i2) {
                        return 0;
                    }
                    int vc = Float.compare(fs[i1],fs[i2]);
                    if (vc != 0) {
                        return - vc;
                    }
                    return i1 > i2 ? 1 : -1;
                }
            });
            BitSet stBS = new BitSet(n);
            //float previous = fs[inds.get(1)];
            //float lastStep = 0;
            for (int sp : inds) {
                stBS.set(sp);
                /*if (previous - fs[sp] < 0) {
						continue;
					}*/
                newBitSets.add((BitSet) stBS.clone());
                //lastStep = previous - fs[sp];
                //previous = fs[sp];
            }
            //System.err.println(this.clusters.getClusterCount());
        }
        return newBitSets;
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
                        is[j] += this.similarityMatrix[k][l];
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
    
    List<BitSet> PhyDstar(SpeciesMapper spm) {
        /* Write the distance matrix to a file,
         *   then use the file as input for PhyD*.
         *   Call PhyD* to construct a tree. */
        File matrix;
        try {
            matrix = new File("tmp.matrix");
            BufferedWriter out = new BufferedWriter(new FileWriter(matrix));

            out.write(spm.getSpeciesCount() + "\n");
            System.out.println(spm.getSpeciesCount() + "\n");
            for (int i = 0; i < spm.getSpeciesCount(); i++) {
                out.write(spm.getSpeciesName(i) + " ");
                System.out.print(spm.getSpeciesName(i) + " ");
                for (int j = 0; j < spm.getSpeciesCount(); j++) {
                    out.write(this.speciesDM[i][j] + " ");
                    System.out.print(this.speciesDM[i][j] + " ");
                }
                System.out.println();
                out.write("\n");
            }
            out.close();
            String[] arg = new String[]{"java","-jar", "PhyDstar.java","-i","tmp.matrix"};
            PhyDstar.main(arg);
        } catch (IOException e) { throw new RuntimeException(); }
        
        /* Write the PhyD* tree back into ASTRAL */
        String newick;
        try {
            File phyDtree = new File("tmp.matrix_bionj.t");
            BufferedReader in = new BufferedReader(new FileReader(phyDtree));
            newick = in.readLine();
            System.out.println(newick);
            in.close();
            matrix.delete();
            phyDtree.delete(); 
        } catch (IOException e) { throw new RuntimeException("Cannot find file: " + e); }
        
        /* Read the newick tree as an actual tree */
        Tree phyDtree = null;
        List<BitSet> ret = new ArrayList<BitSet>();
        try {
            if (newick.length() > 0) {
                newick = newick.replaceAll("\\)[^,);]*", ")");
                NewickReader nr = new NewickReader(new StringReader(newick));
                phyDtree = nr.readTree();
                System.out.println(phyDtree.toNewick());
            }
        } catch (ParseException e) {
            throw new RuntimeException("Failed to Parse Tree: " , e);
        } catch (IOException e) {
            throw new RuntimeException();
        }

        for (TNode node : phyDtree.postTraverse()) {
            if (node.isRoot()) {
                continue;
            } else if (node.isLeaf()) {
                BitSet leaf = new BitSet(n);
                leaf.set(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
                ((STINode)node).setData(leaf);
                System.out.println("Leaf: " + leaf + " " + node.getName());
            } else {
                BitSet newbs = new BitSet(n);
                for (TNode cn : node.getChildren()) {
                    BitSet c = (BitSet) ((STINode)cn).getData();
                    System.out.println("Children: " + c + " " + cn.getName());
                    newbs.or(c);
                }
                ((STINode)node).setData(newbs);
                ret.add(newbs);
            }
        }
        System.out.println(ret);
        return ret;
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
            final float[] is = this.similarityMatrix[i].clone();
            
            sims.add(is);
            range.remove(i);
            TreeSet<Integer> sortColumn = this.sortColumn(range, is);
            range.add(i,i);
            indsBySim.add(sortColumn);
        }

        return upgmaLoop(weights, bsList, indsBySim, sims, n, false);
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
}
