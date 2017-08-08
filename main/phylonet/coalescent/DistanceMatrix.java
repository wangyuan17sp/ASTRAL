package phylonet.coalescent;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.TreeSet;

import phylonet.tree.io.NewickReader;
import phylonet.tree.io.ParseException;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.util.BitSet;

public class DistanceMatrix {
    private float[][] distanceMatrix;
    private List<TreeSet<Integer>> orderedTaxonBySimilarity;
    private Integer n;

    public DistanceMatrix(int n) {
        this.n = n;
    }

    public DistanceMatrix(float[][] from) {
        this.n = from.length;
        this.distanceMatrix = from;
    }

    public int getSize() {
        return n;
    }
 
    public float get(int i, int j) {
        return this.distanceMatrix[i][j];
    }

    List<BitSet> resolveByPhyDstar(List<BitSet> bsList, boolean original) {
        int size = bsList.size();
        float[][] internalMatrix = new float[size][size];
        HashMap<Integer, BitSet> bitMap = new HashMap<Integer, BitSet>(size);

        for (int n = 0; n < 10; n++) {
            int[][] pairLeaves = new int[size][2];
            List<int[]> ones = new ArrayList<int[]>();

            for (int i = 0; i < size; i++) {
                bitMap.put(i, bsList.get(i));
                BitSet bsI = bsList.get(i); 
                int[] bsIone = new int[bsI.cardinality()];
                int c = 0;
                for (int k = bsI.nextSetBit(0); k >= 0; k = bsI.nextSetBit(k + 1)) {
                    bsIone[c++] = k;
                }
                ones.add(bsIone);
                pairLeaves[i][0] = GlobalMaps.random.nextInt(bsI.cardinality());
                pairLeaves[i][1] = GlobalMaps.random.nextInt(bsI.cardinality());
            }

            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    if (i == j) {
                        internalMatrix[i][j] = 0;
                        continue;
                    }
                    
                    float dI1I2 = this.distanceMatrix[ones.get(i)[pairLeaves[i][0]]][ones.get(i)[pairLeaves[i][1]]];
                    float dJ1J2 = this.distanceMatrix[ones.get(j)[pairLeaves[j][0]]][ones.get(j)[pairLeaves[j][1]]];
                    float dI1J1 = this.distanceMatrix[ones.get(i)[pairLeaves[i][0]]][ones.get(j)[pairLeaves[j][0]]];
                    float dI2J2 = this.distanceMatrix[ones.get(i)[pairLeaves[i][1]]][ones.get(j)[pairLeaves[j][1]]];

                    internalMatrix[i][j] += (dI1J1 + dI2J2 - dI1I2 - dJ1J2) / 2;
                }
            }
        }
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                internalMatrix[i][j] /= 10;
            }
        }

        File matrix;
        try {
            matrix = File.createTempFile("internal", ".matrix");
            BufferedWriter out = new BufferedWriter(new FileWriter(matrix));
            out.write(size + "\n");
            for (int i = 0; i < size; i++) {
                out.write(i + " ");
                for (int j = 0; j < size; j++) {
                    out.write(internalMatrix[i][j] + " ");
                }
                out.write("\n");
            }
            out.close();

            String[] arg = new String[]{"java","-jar","PhyDstar.java","-i",matrix.getCanonicalPath()};
            PhyDstar.main(arg);
        } catch (IOException e) { throw new RuntimeException(); }

        Tree phyDtree = generatePhyDstarTree(matrix);

        List<BitSet> ret = new ArrayList<BitSet>();

        /* Extract list of BitSets of leaf-sets */
        for (TNode node : phyDtree.postTraverse()) {
            if (node.isRoot() ) {
                continue;
            } else if (node.isLeaf()) {
                if (original) {
                    BitSet leaf = bitMap.get(Integer.parseInt(node.getName()));
                    ((STINode)node).setData(leaf);
                } else {
                    BitSet leaf = new BitSet(size);
                    leaf.set(Integer.parseInt(node.getName()));
                    ((STINode)node).setData(leaf);
                }
            } else {
                BitSet newbs = new BitSet(n);
                for (TNode cn : node.getChildren()) {
                    BitSet c = (BitSet) ((STINode)cn).getData();
                    newbs.or(c);
                }
                ((STINode)node).setData(newbs);
                ret.add(newbs);
            }
        }
        return ret;
    }

    Tree generatePhyDstarTree(File matrix) {

        /* Write PhyD* tree back into ASTRAL */
        String newick;
        try {
            File phyDtree = new File(matrix.getCanonicalPath() + "_bionj.t");
            BufferedReader in = new BufferedReader(new FileReader(phyDtree));
            newick = in.readLine();
            in.close();
            matrix.delete();
            phyDtree.delete(); 
        } catch (IOException e) { throw new RuntimeException("Cannot find file: " + e); }

        /* Read the newick tree as an actual tree */
        Tree phyDstar_t = null;

        try {
            //            newick = newick.replaceAll("\\)[^,);]*", ")");
            NewickReader nr = new NewickReader(new StringReader(newick));
            phyDstar_t = nr.readTree();

        } catch (ParseException e) {
            throw new RuntimeException("Failed to Parse Tree: " , e);
        } catch (IOException e) {
            throw new RuntimeException();
        }
        return phyDstar_t;
    }

    List<BitSet> PhyDstar(SpeciesMapper spm) {
        /* Write the distance matrix to a file, then use the file as input for PhyD*.
         * Call PhyD* to construct a tree. */
        File matrix;
        try {
            matrix = File.createTempFile("distance", ".matrix");
            BufferedWriter out = new BufferedWriter(new FileWriter(matrix));
            //        out.write(this.speciesDM.length + "\n");
            out.write(n + "\n");
            //        for (int i = 0; i < this.speciesDM.length; i++) {
            for (int i = 0; i < n; i++) {
                out.write(spm.getSpeciesName(i) + " ");
                //            for (int j = 0; j < this.speciesDM.length; j++) {
                for (int j = 0; j < n; j++) {
                    out.write(this.distanceMatrix[i][j] + " ");
                }
                out.write("\n");
            }
            out.close();
            String[] arg = new String[]{"java","-jar","PhyDstar.java","-i",matrix.getCanonicalPath()};
            PhyDstar.main(arg);
        } catch (IOException e) { throw new RuntimeException(); }

        Tree phyDtree = generatePhyDstarTree(matrix);

        List<BitSet> ret = new ArrayList<BitSet>();
        /* Extract list of BitSets of leaf-sets */
        for (TNode node : phyDtree.postTraverse()) {
            if (node.isRoot()) {
                continue;
            } else if (node.isLeaf()) {
                BitSet leaf = new BitSet();
                leaf.set(GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTTaxonIdentifier().taxonId(node.getName()));
                ((STINode)node).setData(leaf);
            } else {
                BitSet newbs = new BitSet(n);
                for (TNode cn : node.getChildren()) {
                    BitSet c = (BitSet) ((STINode)cn).getData();
                    newbs.or(c);
                }
                ((STINode)node).setData(newbs);
                ret.add(newbs);
            }
        }
        return ret;
    }
    
    Iterable<BitSet> getQuadraticBitsets() {
        List<BitSet> newBitSets = new ArrayList<BitSet>();
        ArrayList<Integer> inds = new ArrayList<Integer> (n);
        for (int i = 0; i < n; i++) {
            inds.add(i);
        }
        for (final float[] fs : this.distanceMatrix) {
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
    
    int getBetterSideByFourPoint(int x, int a, int b, int c) {
        double xa = this.distanceMatrix[x][a];
        double xb = this.distanceMatrix[x][b];
        double xc = this.distanceMatrix[x][c];
        double ab = this.distanceMatrix[a][b];
        double ac = this.distanceMatrix[a][c];
        double bc = this.distanceMatrix[b][c];
        double ascore = xa + bc  - (xb + ac); // Note this is similartiy, not distance
        double bscore = xb + ac  - (xa + bc); 
        double cscore = xc + ab - (xb + ac); 
        return ascore >= bscore ?
                ascore >= cscore ? a : c :
                    bscore >= cscore ? b : c;   
    }
    
    private void assureOrderedTaxa () {
        if (this.orderedTaxonBySimilarity == null) {
            this.orderedTaxonBySimilarity = this.sortByDistance(this.distanceMatrix);
        }
    }
    
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

    
    
}
