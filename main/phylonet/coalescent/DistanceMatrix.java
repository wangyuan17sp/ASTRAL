package phylonet.coalescent;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import phylonet.tree.io.NewickReader;
import phylonet.tree.io.ParseException;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.util.BitSet;

public class DistanceMatrix {
    private float[][] distanceMatrix;
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

                internalMatrix[i][j] = (dI1J1 + dI2J2 - dI1I2 - dJ1J2) / 2;
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
                BitSet leaf = new BitSet(n);
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

    /*
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

//        System.out.println(GlobalMaps.taxonIdentifier.taxonCount());
        for (int i = 0; i < n; i++) {
//            System.out.print(GlobalMaps.taxonIdentifier.getTaxonName(i) + " ");
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
//                System.out.print(similarityMatrix[i][j] + " ");         
            }
//            System.out.println();
        }
    } */
}
