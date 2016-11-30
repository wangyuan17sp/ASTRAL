package phylonet.coalescent;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Stack;

import javax.naming.ldap.PagedResultsResponseControl;

import phylonet.coalescent.BipartitionWeightCalculator.Quadrapartition;
import phylonet.coalescent.BipartitionWeightCalculator.Results;
import phylonet.coalescent.BipartitionWeightCalculator.Results2;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;

public class WQInference extends AbstractInference<Tripartition> {
	
	int forceAlg = -1;
	long maxpossible;
	private BipartitionWeightCalculator weightCalculator2;
	
	public WQInference(Options inOptions, List<Tree> trees, List<Tree> extraTrees) {
		super(inOptions, trees, extraTrees);
		
		this.forceAlg = inOptions.getAlg();
	}


	/**
	 * Calculates maximum possible score, to be used for normalization.
	 * @return
	 */
	long calculateMaxPossible() {
		long weight = 0;
		Integer  allsides = null;
		Iterator<STITreeCluster> tit = ((WQDataCollection)this.dataCollection).treeAllClusters.iterator();
		boolean newTree = true;
		
		Deque<Integer> stack = new ArrayDeque<Integer>();
		// TODO: this should not use private stuff from weight calculator. 
		//       redo to use tree objects. 
		for (Integer gtb: ((WQWeightCalculator)this.weightCalculator).geneTreesAsInts()){
			if (newTree) {
				allsides = tit.next().getBitSet().cardinality();
				newTree = false;
			}
			if (gtb >= 0){
				stack.push(1);
			} else if (gtb == Integer.MIN_VALUE) {
				stack.clear();
				newTree = true;
			}  else {
			    ArrayList<Integer> children = new ArrayList<Integer>();
			    Integer newSide = 0;
			    for (int i = gtb; i < 0 ; i++) {
			    	Integer pop = stack.pop();
			        children.add(pop);
			        newSide+=pop;
			    }
			    stack.push(newSide);
                Integer sideRemaining = allsides - newSide;
                if ( sideRemaining !=0) {
                    children.add(sideRemaining);
                }
                for (int i = 0; i < children.size(); i++) {
                	Long a = children.get(i) + 0l;
                    
                    for (int j = i+1; j < children.size(); j++) {
                    	Long b = children.get(j) + 0l;
                        /*if (children.size() > 5) {
                        	if ((side1.s0+side2.s0 == 0? 1 :0) +
                        			(side1.s1+side2.s1 == 0? 1 :0) + 
                        			(side1.s2+side2.s2 == 0? 1:0) > 1)
                        		continue;
                        }
                        */
                        for (int k = j+1; k < children.size(); k++) {
                        	Long c = children.get(k) + 0l;
                            weight += (a+b+c-3) *a*b*c;
                        }
                    }
                }
			}
		}
		return weight/4l;
	}
	
	void initializeWeightCalculator() {
		((WQWeightCalculator)this.weightCalculator).setupGeneTrees(this);
		if (this.forceAlg == 2) {
			((WQWeightCalculator)this.weightCalculator).useSetWeightsAlgorithm();
		} 

		this.weightCalculator.initializeWeightContainer(
				this.trees.size() *  GlobalMaps.taxonIdentifier.taxonCount() * 2);
	}
	
	/**
	 * Score first computes the quartet scores and the calls
	 * scoreBranches to annotate branches (if needed). 
	 */
	public double scoreSpeciesTreeWithGTLabels(Tree st, boolean initialize) {

		if (initialize) {
			mapNames();
	
			IClusterCollection clusters = newClusterCollection();
	
	
			this.dataCollection = newCounter(clusters);
			weightCalculator = newWeightCalculator();
	
			WQDataCollection wqDataCollection = (WQDataCollection) this.dataCollection;
			wqDataCollection.preProcess(this);
			this.initializeWeightCalculator();			
			this.maxpossible = this.calculateMaxPossible();
			System.err.println("Number of quartet trees in the gene trees: "+this.maxpossible);
	
			//System.err.println(this.maxpossible);
		}
		
		Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
		long sum = 0l;

		for (TNode node: st.postTraverse()) {
			if (node.isLeaf()) {
				String nodeName = node.getName(); //GlobalMaps.TaxonNameMap.getSpeciesName(node.getName());

				STITreeCluster cluster = new STITreeCluster();
				Integer taxonID = GlobalMaps.taxonIdentifier.taxonId(nodeName);
				cluster.addLeaf(taxonID);

				stack.add(cluster);

			} else {
				ArrayList<STITreeCluster> childbslist = new ArrayList<STITreeCluster>();
				BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
				for (TNode child: node.getChildren()) {
					STITreeCluster pop = stack.pop();
					childbslist.add(pop);
					bs.or(pop.getBitSet());
				}
				
				STITreeCluster cluster = new STITreeCluster();
				cluster.setCluster((BitSet) bs.clone());

				//((STINode)node).setData(new GeneTreeBitset(node.isRoot()? -2: -1));
				stack.add(cluster);


				STITreeCluster remaining = cluster.complementaryCluster();
				if (remaining.getClusterSize() != 0) {
					childbslist.add(remaining);
				}
				if (childbslist.size() > 3) {
					for (STITreeCluster chid :childbslist) {
						System.err.print(chid.getClusterSize()+" ");
					}
					System.err.println(" (polytomy)");
					if (this.getBranchAnnotation() % 2 == 0) {
						continue;
					}
				}
				
				for (int i = 0; i < childbslist.size(); i++) {
					for (int j = i+1; j < childbslist.size(); j++) {
						for (int k = j+1; k < childbslist.size(); k++) {
							Tripartition trip = new Tripartition(childbslist.get(i),  childbslist.get(j), childbslist.get(k));
							Long s = weightCalculator.getWeight(trip, null);
							sum += s;
						}
					}					       
				}
			}
		}
		

		System.err.println("Quartet score is: " + sum/4l);
		System.err.println("Normalized quartet score is: "+ (sum/4l+0.)/this.maxpossible);
		//System.out.println(st.toNewickWD());

		if (this.getBranchAnnotation() == 0){
			for (TNode n: st.postTraverse()) {
				((STINode) n).setData(null);
			}
		} else {
			if (this.getDepth() == 0){
				double logscore = this.scoreBranches(st);
				if (this.getBranchAnnotation() % 12 == 0) {
					System.err.println("log local posterior: "+logscore);
					return logscore;
				}
			}
			else {
				scoreBranches(st, this.getDepth());
			}
		}
		return (sum/4l+0.)/this.maxpossible;
		
	}
	
	
	private boolean skipNode (TNode node) {
		return 	node.isLeaf() || node.isRoot() || node.getChildCount() > 2 || 
				(node.getParent().getChildCount() >3) ||
				(node.getParent().getChildCount() >2 && !node.getParent().isRoot())||
				((node.getParent().isRoot() && node.getParent().getChildCount() == 2));
	}
	
	private class NodeData {
		double mainfreq, alt1freqs, alt2freqs;
		long quartcount;
		double effn ;
		Quadrapartition [] quads  = {null, null, null};
		STBipartition[] bipartitions = {null, null, null};
		double post;
		Posterior postQ1 = null;
		Posterior postQ2 = null;
		Posterior postQ3 = null;
		String data;
		void setData(String s){
			this.data = s;
		}
		void setString(int i){
			Double f1 = mainfreq;
			Double f2 = alt1freqs;
			Double f3 = alt2freqs;
			Double effni = effn;
			
			if (i == 0){
				this.setData(null);
			} else if (i == 1){
				this.setData(df.format((f1+.0)/effni*100));
			} else if (i == 10) {
				df.setMaximumFractionDigits(5);
				this.setData(df.format(this.postQ1.getPost()));
			} else {
				double pQ1 = this.postQ1.getPost();
				
				if (i == 3 || i == 12) {
					this.setData(df.format(postQ1.getPost()));
				} else if (i % 2 == 0) {
					double pQ2 = postQ2.getPost();
					double pQ3 = postQ3.getPost();
					if (i == 2)
						this.setData("'[q1="+(f1)/effni+";q2="+(f2)/effni+";q3="+(f3)/effni+
								 ";f1="+f1+";f2="+f2+";f3="+f3+
								 ";pp1="+pQ1+";pp2="+pQ2+";pp3="+pQ3+
								 ";QC="+this.quartcount+";EN="+effni+"]'");
					else if (i == 4) {
						this.setData("'[pp1="+df.format(pQ1)+";pp2="+df.format(pQ2)+";pp3="+df.format(pQ3)+"]'");
					} else if (i == 6){
						this.setData(df.format(pQ1));
						this.setData(this.quads[0] +
								" [" + this.bipartitions[0].toString2() + "] : " + pQ1 + " ** f1 = " + f1 +
								" f2 = " + f2 + " f3 = " + f3 + " EN = " + effni + " **\n" +  this.quads[1] +
								" [" + this.bipartitions[1].toString2() + "] : " + pQ2 + " ** f1 = " + f2 +
								" f2 = " + f1 + " f3 = " + f3 + " EN = " + effni + " **\n"+this.quads[2] +
								" [" + this.bipartitions[2].toString2() + "] : " + pQ3 + " ** f1 = " + f3 +
								" f2 = " + f1 + " f3 = " + f2 + " EN = " + effni + " **\n");
					}  else if (i == 8){
						this.setData(
								"'[q1="+df.format((f1)/effni)+
								 ";q2="+df.format((f2)/effni)+
								 ";q3="+df.format((f3)/effni)+"]'");
					}
				}
			}
		}
		String toString2(int i){
			Double f1 = mainfreq;
			Double f2 = alt1freqs;
			Double f3 = alt2freqs;
			Double effni = effn;
			
			if (i == 0){
				return null;
			} else if (i == 1){
				return df.format((f1+.0)/effni*100);
			} else if (i == 10) {
				df.setMaximumFractionDigits(5);
				return df.format(this.postQ1.getPost());
			} else {
				double pQ1 = this.postQ1.getPost();
				
				if (i == 3 || i == 12) {
					return df.format(pQ1);
				} else if (i % 2 == 0) {
					double pQ2 = postQ2.getPost();
					double pQ3 = postQ3.getPost();
					if (i == 2)
						return "'[q1="+(f1)/effni+";q2="+(f2)/effni+";q3="+(f3)/effni+
								 ";f1="+f1+";f2="+f2+";f3="+f3+
								 ";pp1="+pQ1+";pp2="+pQ2+";pp3="+pQ3+
								 ";QC="+this.quartcount+";EN="+effni+"]'";
					else if (i == 4) {
						return "'[pp1="+df.format(pQ1)+";pp2="+df.format(pQ2)+";pp3="+df.format(pQ3)+"]'";
					} else if (i == 6){
						return df.format(pQ1);
					}  else if (i == 8){
						return "'[q1="+df.format((f1)/effni)+
								 ";q2="+df.format((f2)/effni)+
								 ";q3="+df.format((f3)/effni)+"]'";
					}
				}
			}
			return null;
		}
		
		String toString2(){
			return this.data;
		}
		
	}
	private ArrayList<STITreeCluster> listClustersBelowNode(TNode node, int depth){
		ArrayList<STITreeCluster> retClusters = new ArrayList<STITreeCluster>();
		STINode n = (STINode) node;
		if (depth == 0 || node.isLeaf()){
				retClusters.add((STITreeCluster) n.getData());
				return retClusters;
		}
		for(TNode child: node.getChildren()){
			retClusters.addAll(listClustersBelowNode(child, depth - 1));
		}
		return retClusters;
		
	}	
	
	private ArrayList<STITreeCluster> listClustersAboveParentNode(TNode node, int maxDist){
		ArrayList<STITreeCluster> retClusters = new ArrayList<STITreeCluster>();
		STINode n = (STINode) node;
		if (maxDist == 0 || node.isLeaf()) {
			STITreeCluster tmp = ((STITreeCluster) n.getData()).complementaryCluster();
			retClusters.add(tmp);
			return retClusters;
		}
		Iterator<STINode> pcit = n.getParent().getChildren().iterator();
		TNode child = (TNode) pcit.next();
		if (child == node) child = pcit.next();
		retClusters.addAll(listClustersBelowNode(child, maxDist - 1));
		if (n.getParent().isRoot()) {
			return retClusters;
		}
		else {
			retClusters.addAll(listClustersAboveParentNode(n.getParent(),maxDist - 1));
		}
		return retClusters;
	}
	
	private boolean toAnnotateAltTopologies() {
		if (this.getBranchAnnotation() != 0 && this.getBranchAnnotation() != 1 && this.getBranchAnnotation() != 10 && 
				this.getBranchAnnotation() != 3 && this.getBranchAnnotation() != 12 && this.getBranchAnnotation() % 2 == 0) {
			return true;
		}
		return false;
	}
	
	private void setLocalPP(NodeData nd){
		double f1  = nd.mainfreq;
		double f2  = nd.alt1freqs;
		double f3  = nd.alt2freqs;
		long quarc = nd.quartcount;
		double effni = nd.effn + 0.0;
		
		if ( Math.abs((f1+f2+f3) - effni) > 0.01 ) {
			effni = f1 + f2 + f3;
		}
		
		
		Posterior post = new Posterior(
				f1,f2,f3,(double)effni, options.getLambda());
		System.err.println("the computed localpp is:" + post + " " + f1 + " " + f2 + " " + f3 + " " + effni);
		nd.postQ1 = post;
		nd.post = post.getPost();
		
		if (toAnnotateAltTopologies()) {
			post = new Posterior(f2,f1,f3,(double)effni, options.getLambda());
			nd.postQ2 = post;
			post =  new Posterior(f3,f1,f2,(double)effni, options.getLambda());
			nd.postQ3 = post;
		}
	}
	
	private void scoreBranches(Tree st, int depth){
		System.out.println(depth);
		System.out.println(options.getBranchannotation());
		weightCalculator = new BipartitionWeightCalculator(this,((WQWeightCalculator)this.weightCalculator).geneTreesAsInts());
		weightCalculator2 = (BipartitionWeightCalculator) weightCalculator;
				
		/**
		 * Add bitsets to each node for all taxa under it. 
		 * Bitsets are saved in nodes "data" field
		 */
		Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;
			if (node.isLeaf()) {
				String nodeName = node.getName(); //GlobalMaps.TaxonNameMap.getSpeciesName(node.getName());

				STITreeCluster cluster = new STITreeCluster();
				Integer taxonID = GlobalMaps.taxonIdentifier.taxonId(nodeName);
				cluster.addLeaf(taxonID);

				stack.add(cluster);
				node.setData(cluster);

			} else {
				ArrayList<STITreeCluster> childbslist = new ArrayList<STITreeCluster>();
				BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
				for (TNode child: n.getChildren()) {
					STITreeCluster pop = stack.pop();
					childbslist.add(pop);
					bs.or(pop.getBitSet());
				}
				
				STITreeCluster cluster = new STITreeCluster();
				cluster.setCluster((BitSet) bs.clone());

				//((STINode)node).setData(new GeneTreeBitset(node.isRoot()? -2: -1));
				stack.add(cluster);
				node.setData(cluster);
			}
		}
		
		stack = new Stack<STITreeCluster>();
		Map<TNode,NodeData> finalNodeDataList = new HashMap<TNode,NodeData>();
		for (TNode n: st.postTraverse()) {
			if (n.isLeaf() || n.getParent() == null || n.getParent().getParent() == null
					|| n.getChildCount() > 2) {
				continue;
			} else {
				finalNodeDataList.put(n,scoreSTNode(depth, n));
			}
		}
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;

			if (finalNodeDataList.containsKey(n)) {
				node.setData(finalNodeDataList.get(node).toString2(this.getBranchAnnotation()));
			} else {
				node.setData(null);
			}
			
		}
	}


	private NodeData scoreSTNode(int depth, TNode n) {
		
		ArrayList<NodeData> nodeDataList = new ArrayList<NodeData>();
		STINode node = (STINode) n;

		/**
		 * 1. Create quadripartion
		 */
		ArrayList<STITreeCluster>[] cllists = clustersAtDepth(depth, node);
		ArrayList<STITreeCluster> belowClusters1 = cllists[0];
		ArrayList<STITreeCluster> belowClusters2 = cllists[1];
		ArrayList<STITreeCluster> sisterClusters = cllists[2];
		ArrayList<STITreeCluster> remainingClusters = cllists[3];

		STITreeCluster[] topcls = summerizeClusters(belowClusters1, belowClusters2, sisterClusters, remainingClusters);
		
		STITreeCluster belowCl1 = topcls[0];
		STITreeCluster belowCl2 = topcls[1];
		STITreeCluster sisterCl = topcls[2];
		STITreeCluster remainingCl = topcls[3];

		
		System.err.println(belowClusters1.size() +" "+belowClusters2.size() + " "+
				sisterClusters.size() + " "+ remainingClusters.size());

		if (belowClusters1.size() +belowClusters2.size() + 
				sisterClusters.size() +  remainingClusters.size() > 2<<(depth)) {
			System.err.println("polytomy; seeing more than " + (2<<(depth))+ " "+depth);
			return null;
		}
		
		NodeData summerizedNd = new NodeData();
		if (this.options.getSumMethod() == 0 || this.options.getSumMethod() == 3) {
				nodeDataList = ScoreAllPossibleQuadripartitions(belowClusters1, belowClusters2, sisterClusters, remainingClusters);
				if (nodeDataList.size()==0) {
					throw new RuntimeException("Hmm, this shouldn't happen; "+nodeDataList);
				}
				if (this.options.getSumMethod() == 0) {
					summerizedNd = findNDWithMinimumLocalPP(nodeDataList);
				}
				if (this.options.getSumMethod() == 3) {
					summerizedNd = findNDWithAverageLoalPP(nodeDataList, belowCl1, belowCl2, sisterCl, remainingCl);
				}
		}
		else if (this.options.getSumMethod() == 1 || this.options.getSumMethod() == 4) {
			nodeDataList = ScoreOFESQuadripartitions(belowClusters1, belowCl1, belowClusters2, belowCl2, sisterClusters, sisterCl,  remainingClusters, remainingCl);
			if (nodeDataList.size()==0) {
				throw new RuntimeException("Hmm, this shouldn't happen; "+nodeDataList);
			}
			if (this.options.getSumMethod() == 1) {
				summerizedNd = findNDWithMinimumLocalPP(nodeDataList);
			}
			if (this.options.getSumMethod() == 4) {
				summerizedNd = findNDWithAverageLoalPP(nodeDataList, belowCl1, belowCl2, sisterCl, remainingCl);
			}
		}
		else if (this.options.getSumMethod() == 2 || this.options.getSumMethod() == 5) { 	
			nodeDataList = ScoreAllPlusOFESQuadripartitions(belowClusters1, belowCl1, belowClusters2, belowCl2, sisterClusters, sisterCl,  remainingClusters, remainingCl);
			if (nodeDataList.size()==0) {
				throw new RuntimeException("Hmm, this shouldn't happen; "+nodeDataList);
			}
			if (this.options.getSumMethod() == 2) {
				summerizedNd = findNDWithMinimumLocalPP(nodeDataList);
			}
			if (this.options.getSumMethod() == 5) {
				summerizedNd = findNDWithAverageLoalPP(nodeDataList, belowCl1, belowCl2, sisterCl, remainingCl);
			}
		}
		else {
			throw new RuntimeException("Hmm, bad summerization method option; "+this.options.getSumMethod());
		}
		
		
		STITreeCluster cluster = new STITreeCluster(belowCl1);
		cluster.getBitSet().or(belowCl2.getBitSet());
		STBipartition bmain = new STBipartition(cluster , cluster.complementaryCluster());
		
		STITreeCluster c1plussis = new STITreeCluster(belowCl1);
		c1plussis.getBitSet().or(sisterCl.getBitSet());
		STBipartition b2 = new STBipartition(c1plussis, c1plussis.complementaryCluster());
		
		STITreeCluster c1plusrem = new STITreeCluster(belowCl1);
		c1plusrem.getBitSet().or(remainingCl.getBitSet());
		STBipartition b3 = new STBipartition(c1plusrem, c1plusrem.complementaryCluster());	
		
		STBipartition[] biparts = new STBipartition[] {bmain, b2, b3};	
		summerizedNd.bipartitions = biparts;

		reorderClusters(summerizedNd, bmain, b2, b3);
		
		
		double bl = summerizedNd.postQ1.branchLength();
		
		node.setParentDistance(bl);
		summerizedNd.setString(this.getBranchAnnotation());		
		System.err.println(summerizedNd.data);
		
		return summerizedNd;

	}
	
	
	

	private void reorderClusters(NodeData node, STBipartition bmain,
			STBipartition b2, STBipartition b3) {
		// TODO Auto-generated method stub
		reorderACluster(node.quads[0], bmain);
		if (toAnnotateAltTopologies()) {
			reorderACluster(node.quads[1], b2);
			reorderACluster(node.quads[2], b3);
		}
	}


	private void reorderACluster(Quadrapartition quad,
			STBipartition bipart) {
		// TODO Auto-generated method stub
		STITreeCluster tmpcluster12 = new STITreeCluster(quad.cluster1);
		tmpcluster12.getBitSet().or(quad.cluster2.getBitSet());
		
		STITreeCluster tmpcluster34 = new STITreeCluster(quad.cluster3);
		tmpcluster34.getBitSet().or(quad.cluster4.getBitSet());

		
		if (!bipart.cluster1.getBitSet().contains(tmpcluster12.getBitSet()) || 
				!bipart.cluster2.getBitSet().contains(tmpcluster34.getBitSet())) {
			swapClusters12AndClusters34(quad);
			checkIfQuadMatchesBipart(quad, bipart);
		}
		
	}


	private void checkIfQuadMatchesBipart(Quadrapartition quad , STBipartition bipart) {
		// TODO Auto-generated method stub
		STITreeCluster tmpcluster12 = new STITreeCluster(quad.cluster1);
		STITreeCluster tmpcluster34 = new STITreeCluster(quad.cluster3);
		
		tmpcluster12.getBitSet().or(quad.cluster2.getBitSet());
		tmpcluster34.getBitSet().or(quad.cluster4.getBitSet());
		
		if (!bipart.cluster1.getBitSet().contains(tmpcluster12.getBitSet()) || 
				!bipart.cluster2.getBitSet().contains(tmpcluster34.getBitSet())) {

			throw new RuntimeException("Hmm, this shouldn't happen; "+ tmpcluster12.getBitSet().nextSetBit(0) + " " + tmpcluster34.getBitSet().nextSetBit(0));
		}
	}


	private void swapClusters12AndClusters34(Quadrapartition quad) {
		// TODO Auto-generated method stub
		STITreeCluster tmp1 = new STITreeCluster(quad.cluster3);
		STITreeCluster tmp2 = new STITreeCluster(quad.cluster4);
		quad.cluster3 = quad.cluster1;
		quad.cluster4 = quad.cluster2;
		quad.cluster1 = tmp1;
		quad.cluster2 = tmp2;
	}
	


	private void extracted(NodeData node) {
		STITreeCluster tmp1 = new STITreeCluster(node.quads[1].cluster3);
		STITreeCluster tmp2 = new STITreeCluster(node.quads[1].cluster4);
		node.quads[1].cluster3 = node.quads[1].cluster1;
		node.quads[1].cluster4 = node.quads[1].cluster2;
		node.quads[1].cluster1 = tmp1;
		node.quads[1].cluster2 = tmp2;
	}


	private STITreeCluster[] summerizeClusters(
			ArrayList<STITreeCluster> belowClusters1,
			ArrayList<STITreeCluster> belowClusters2,
			ArrayList<STITreeCluster> sisterClusters,
			ArrayList<STITreeCluster> remainingClusters) {
		
		STITreeCluster cluster1  = new STITreeCluster();
		STITreeCluster cluster2  = new STITreeCluster();
		STITreeCluster sister    = new STITreeCluster();
		STITreeCluster remaining = new STITreeCluster();
		
		for (STITreeCluster c: belowClusters1) {
			cluster1.getBitSet().or(c.getBitSet());
		}
		
		for (STITreeCluster c: belowClusters2) {
			cluster2.getBitSet().or(c.getBitSet());
		}
		
		for (STITreeCluster c: sisterClusters) {
			sister.getBitSet().or(c.getBitSet());
		}
		
		for (STITreeCluster c: remainingClusters) {
			remaining.getBitSet().or(c.getBitSet());
		}
		
		return new STITreeCluster [] {cluster1, cluster2, sister, remaining};
	}


	private NodeData findNDWithAverageLoalPP(ArrayList<NodeData> nodeDataList, STITreeCluster c1, STITreeCluster c2, STITreeCluster sister, STITreeCluster remaining) {
		NodeData criticalNd= new NodeData();
		for(NodeData ndI: nodeDataList){
				if (ndI == null) {
					throw new RuntimeException("Hmm, this shouldn't happen; "+ndI);
				}

				criticalNd.mainfreq += ndI.mainfreq;

				criticalNd.alt1freqs += ndI.alt1freqs;
				criticalNd.alt2freqs += ndI.alt2freqs;				
				criticalNd.effn += ndI.effn;
				criticalNd.quartcount += ndI.quartcount;

		} 
		criticalNd.mainfreq /= nodeDataList.size();
		criticalNd.alt1freqs /= nodeDataList.size();
		criticalNd.alt2freqs /= nodeDataList.size();
		criticalNd.effn /= nodeDataList.size();
		criticalNd.quartcount /= nodeDataList.size();

		criticalNd.postQ1 = new Posterior(criticalNd.mainfreq, 
				criticalNd.alt1freqs, criticalNd.alt2freqs, 
				criticalNd.effn, options.getLambda());
		
		criticalNd.postQ2 = new Posterior(criticalNd.alt1freqs, 
				criticalNd.mainfreq, criticalNd.alt2freqs, 
				criticalNd.effn, options.getLambda());
		
		criticalNd.postQ3 = new Posterior(criticalNd.alt2freqs, 
				criticalNd.mainfreq, criticalNd.alt1freqs, 
				criticalNd.effn, options.getLambda());
		criticalNd.quads[0] = weightCalculator2.new Quadrapartition
				(c1, c2, sister, remaining, true);
		
		criticalNd.quads[1] = weightCalculator2.new Quadrapartition
				(c1, sister, c2, remaining, true);
		criticalNd.quads[2] = weightCalculator2.new Quadrapartition
				(c1, remaining, c2, sister, true);
		return criticalNd;
	}


	private NodeData findNDWithMinimumLocalPP(ArrayList<NodeData> nodeDataList) {
		double minPostQ1 = Double.MAX_VALUE;
		double minPostQ2 = Double.MAX_VALUE;
		double minPostQ3 = Double.MAX_VALUE;
		NodeData criticalNd= new NodeData();
		for(NodeData ndI: nodeDataList){
				if (ndI == null) {
					throw new RuntimeException("Hmm, this shouldn't happen; "+ndI);
				}
				if (ndI.postQ1.getPost() < minPostQ1) {
					minPostQ1 = ndI.postQ1.getPost();
					criticalNd.postQ1 = ndI.postQ1;
					criticalNd.quads[0] = ndI.quads[0];
					criticalNd.mainfreq = ndI.mainfreq;
					criticalNd.effn = ndI.effn;
					criticalNd.quartcount = ndI.quartcount;
				}
				if (toAnnotateAltTopologies()) {
					if (ndI.postQ2.getPost() < minPostQ2) {
						minPostQ2 = ndI.postQ2.getPost();
						criticalNd.postQ2 = ndI.postQ2;
						criticalNd.quads[1] = ndI.quads[1];
						criticalNd.alt1freqs = ndI.alt1freqs;
					}
				
					if (ndI.postQ3.getPost() < minPostQ3) {
						minPostQ3 = ndI.postQ3.getPost();
						criticalNd.postQ3 = ndI.postQ3;
						criticalNd.quads[2] = ndI.quads[2];
						criticalNd.alt2freqs = ndI.alt2freqs;
					}
				}

		} 

		return criticalNd;
	}


	private ArrayList<NodeData> ScoreOFESQuadripartitions(
			ArrayList<STITreeCluster> belowClusters1, STITreeCluster belowCl1,
			ArrayList<STITreeCluster> belowClusters2, STITreeCluster belowCl2,
			ArrayList<STITreeCluster> sisterClusters, STITreeCluster sisterCl,
			ArrayList<STITreeCluster> remainingClusters,
			STITreeCluster remainingCl) {
		ArrayList<NodeData> nodeDataList = new ArrayList<NodeData>();
		
			
		
		for (STITreeCluster c1: belowClusters1) {
	
			nodeDataList.add(ScoreAQuadripartion(c1, belowCl2, sisterCl, remainingCl));
		}
		for (STITreeCluster c2: belowClusters2) {
			
			nodeDataList.add(ScoreAQuadripartion(belowCl1, c2, sisterCl, remainingCl));
		}
		for (STITreeCluster sister: sisterClusters) {
			
			nodeDataList.add(ScoreAQuadripartion(belowCl1, belowCl2, sister, remainingCl));
		}
		for (STITreeCluster remaining: remainingClusters) {
			
			nodeDataList.add(ScoreAQuadripartion(belowCl1, belowCl2, sisterCl, remaining));
		}
			
		return nodeDataList;
	}


	private ArrayList<NodeData> ScoreAllPlusOFESQuadripartitions(
			ArrayList<STITreeCluster> belowClusters1, STITreeCluster belowCl1,
			ArrayList<STITreeCluster> belowClusters2, STITreeCluster belowCl2,
			ArrayList<STITreeCluster> sisterClusters, STITreeCluster sisterCl,
			ArrayList<STITreeCluster> remainingClusters,
			STITreeCluster remainingCl) {
		ArrayList<NodeData> nodeDataList = new ArrayList<NodeData>();
		for (STITreeCluster c1: belowClusters1) {
	
			nodeDataList.add(ScoreAQuadripartion(c1, belowCl2, sisterCl, remainingCl));
		}
		for (STITreeCluster c2: belowClusters2) {
			
			nodeDataList.add(ScoreAQuadripartion(belowCl1, c2, sisterCl, remainingCl));
		}
		for (STITreeCluster sister: sisterClusters) {
			
			nodeDataList.add(ScoreAQuadripartion(belowCl1, belowCl2, sister, remainingCl));
		}
		for (STITreeCluster remaining: remainingClusters) {
			
			nodeDataList.add(ScoreAQuadripartion(belowCl1, belowCl2, sisterCl, remaining));
		}
			
		for (STITreeCluster c1: belowClusters1){
			for (STITreeCluster c2: belowClusters2){
				for (STITreeCluster sister: sisterClusters){					
					for (STITreeCluster remaining: remainingClusters){
						nodeDataList.add(ScoreAQuadripartion(c1, c2, sister, remaining));
					}
				}
			}
		}
		
		return nodeDataList;
		
	}



	private ArrayList<STITreeCluster>[] clustersAtDepth(int depth, STINode node) {

		ArrayList<STITreeCluster> belowClusters1;
		ArrayList<STITreeCluster> belowClusters2;
		ArrayList<STITreeCluster> sisterClusters;
		ArrayList<STITreeCluster> remainingClusters;
			
		Iterator<STINode> pcit = node.getChildren().iterator();
		TNode pc = (TNode) pcit.next();
		belowClusters1 = listClustersBelowNode(pc, depth-1);
		pc = pcit.next();
		belowClusters2 = listClustersBelowNode(pc, depth-1);
		pcit = node.getParent().getChildren().iterator();
		pc = pcit.next();
		while (pc == node) pc = pcit.next();
		sisterClusters = listClustersBelowNode(pc, depth-1);
		remainingClusters = listClustersAboveParentNode(node.getParent(),depth-1);
			
		return new ArrayList [] {belowClusters1,belowClusters2, sisterClusters, remainingClusters};
	}



	private ArrayList<NodeData> ScoreAllPossibleQuadripartitions(
			ArrayList<STITreeCluster> belowClusters1,
			ArrayList<STITreeCluster> belowClusters2,
			ArrayList<STITreeCluster> sisterClusters,
			ArrayList<STITreeCluster> remainingClusters) {
		ArrayList<NodeData> nodeDataList = new ArrayList<NodeData>();
		for (STITreeCluster c1: belowClusters1){
			for (STITreeCluster c2: belowClusters2){
				for (STITreeCluster sister: sisterClusters){					
					for (STITreeCluster remaining: remainingClusters){
						nodeDataList.add(ScoreAQuadripartion(c1, c2, sister, remaining));
					}
				}
			}
		}
		return nodeDataList;
	}


	private NodeData ScoreAQuadripartion(STITreeCluster c1, STITreeCluster c2,
			STITreeCluster sister, STITreeCluster remaining) {
		
		//TODO: figure out this
		if (this.getBranchAnnotation() == 7){
			if (remaining.getClusterSize() != 0 && sister.getClusterSize() != 0 && c2.getClusterSize() != 0 && c1.getClusterSize() != 0 ){
				System.err.print(c1.toString()+c2.toString()+"|"+sister.toString()+remaining.toString()+"\n");
			}
		}
		NodeData nd = new NodeData();

		
		Quadrapartition quad = weightCalculator2.new Quadrapartition
				(c1,  c2, sister, remaining, false);
		
		
		Results s = weightCalculator2.getWeight(quad);
		nd.mainfreq = s.qs[0];
		nd.alt1freqs = s.qs[2];	
		nd.alt2freqs = s.qs[1];
		nd.effn = (double) s.effn + 0.0;		
		nd.quartcount = (c1.getClusterSize()+0l)
				* (c2.getClusterSize()+0l)
				* (sister.getClusterSize()+0l)
				* (remaining.getClusterSize()+0l);
		
		
		nd.quads[0] = weightCalculator2.new Quadrapartition
				(c1, c2, sister, remaining, true);		
		nd.quads[1] = weightCalculator2.new Quadrapartition
				(c1, sister, c2, remaining, true);
		nd.quads[2] = weightCalculator2.new Quadrapartition
				(c1, remaining, c2, sister, true);
		
		this.setLocalPP(nd);
		return nd;
	}
	

	
	
	/**
	 * Annotates the species tree branches with support, branch length, etc. 
	 * @param st
	 * @return
	 */
	private double scoreBranches(Tree st) {

		double ret = 0;
		
		weightCalculator = new BipartitionWeightCalculator(this,((WQWeightCalculator)this.weightCalculator).geneTreesAsInts());
		
		BipartitionWeightCalculator weightCalculator2 = (BipartitionWeightCalculator) weightCalculator;
		WQDataCollection wqDataCollection = (WQDataCollection) this.dataCollection;
		//wqDataCollection.initializeWeightCalculator(this);
		
		/**
		 * Add bitsets to each node for all taxa under it. 
		 * Bitsets are saved in nodes "data" field
		 */
		Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;
			if (node.isLeaf()) {
				String nodeName = node.getName(); //GlobalMaps.TaxonNameMap.getSpeciesName(node.getName());

				STITreeCluster cluster = new STITreeCluster();
				Integer taxonID = GlobalMaps.taxonIdentifier.taxonId(nodeName);
				cluster.addLeaf(taxonID);

				stack.add(cluster);
				node.setData(cluster);

			} else {
				ArrayList<STITreeCluster> childbslist = new ArrayList<STITreeCluster>();
				BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
				for (TNode child: n.getChildren()) {
					STITreeCluster pop = stack.pop();
					childbslist.add(pop);
					bs.or(pop.getBitSet());
				}
				
				STITreeCluster cluster = new STITreeCluster();
				cluster.setCluster((BitSet) bs.clone());

				//((STINode)node).setData(new GeneTreeBitset(node.isRoot()? -2: -1));
				stack.add(cluster);
				node.setData(cluster);
			}
		}
		stack = new Stack<STITreeCluster>();
		
		
		/**
		 * For each node,
		 *   1. create three quadripartitoins for the edge above it
		 *   2. score the quadripartition
		 *   3. save the scores in a list for annotations in the next loop
		 */
		Map<TNode,NodeData> nodeDataList = new HashMap<TNode,NodeData>();
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;
			if (node.isLeaf()) {
				stack.push((STITreeCluster) node.getData());
			} else {
				/**
				 * 1. Create quadripartion
				 */
				NodeData nd = new NodeData();
				
				
				STITreeCluster cluster = (STITreeCluster) node.getData();				
				if (skipNode(node) ) {
					for (int i =0; i< node.getChildCount(); i++) {
						stack.pop();
					}
					stack.push(cluster);
					continue;
				}
				
				
				STITreeCluster c2 = stack.pop();
				STITreeCluster c1 = stack.pop();
				stack.push(cluster);
				
				STITreeCluster sister;
				STITreeCluster remaining;
				Iterator<STINode> pcit = node.getParent().getChildren().iterator();
				STINode pc = pcit.next();
				if ( pc == n ) pc = pcit.next(); 
				sister = (STITreeCluster)pc.getData();
				if (node.getParent().isRoot() && node.getParent().getChildCount() == 3) {
					pc = pcit.next();
					if (pc == n) pc = pcit.next(); 
					remaining = (STITreeCluster)pc.getData();;					
				} /* else if (node.getParent().isRoot() && node.getParent().getChildCount() == 2) {
					continue;
				} */ 
				else {
					remaining = ((STITreeCluster)node.getParent().getData()).complementaryCluster();
				}
				
				if (this.getBranchAnnotation() == 7){
					if (remaining.getClusterSize() != 0 && sister.getClusterSize() != 0 && c2.getClusterSize() != 0 && c1.getClusterSize() != 0 ){
						System.err.print(c1.toString()+c2.toString()+"|"+sister.toString()+remaining.toString()+"\n");
					}
				}
				
				/**
				 * 2. Scores all three quadripartitoins
				 */
				
				
				
				if (nd.effn < 20) {
					if (!GlobalMaps.taxonNameMap.getSpeciesIdMapper().isSingleSP(cluster.getBitSet()))
						System.err.println("You may want to ignore posterior probabilities and other statistics related to the following "
								+ "branch branch because the effective number of genes impacting it is only "+ nd.effn +
							":\n\t" +
							GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTClusterForGeneCluster(cluster));
				}
				
				Quadrapartition quadmain = weightCalculator2.new Quadrapartition
						(c1,  c2, sister, remaining, true);
				
				Quadrapartition quad2 = weightCalculator2.new Quadrapartition
						(c1, sister, c2, remaining, true);
				
				Quadrapartition quad3 = weightCalculator2.new Quadrapartition
						(c1, remaining, c2, sister, true);
				
				
				Quadrapartition[] threequads = new Quadrapartition [] {quadmain, quad2, quad3};
				
				Quadrapartition quad = weightCalculator2.new Quadrapartition
						(c1,  c2, sister, remaining, false);
				
				Results s = weightCalculator2.getWeight(quad);
				nd.mainfreq = s.qs[0];
				nd.effn = s.effn;
				nd.alt1freqs = s.qs[2];				
				nd.alt2freqs = s.qs[1];
				nd.quartcount= (c1.getClusterSize()+0l)
						* (c2.getClusterSize()+0l)
						* (sister.getClusterSize()+0l)
						* (remaining.getClusterSize()+0l);
			
				STITreeCluster c1plussis = new STITreeCluster();
				c1plussis.setCluster((BitSet) c1.getBitSet().clone());
				
				c1plussis.getBitSet().or(sister.getBitSet());
				STITreeCluster c1plusrem = new STITreeCluster();
				
				c1plusrem.setCluster((BitSet) c1.getBitSet().clone());
				c1plusrem.getBitSet().or(remaining.getBitSet());
					
				STBipartition bmain = new STBipartition(cluster, cluster.complementaryCluster());
				STBipartition b2 = new STBipartition(c1plussis, c1plussis.complementaryCluster());
				STBipartition b3 = new STBipartition(c1plusrem, c1plusrem.complementaryCluster());
	
				STBipartition[] biparts = new STBipartition[] {bmain, b2, b3};
				nd.quads = threequads;
				nd.bipartitions = biparts;

				reorderClusters(nd, bmain, b2, b3);
				nodeDataList.put(node, nd);
			}
		}
		
		/**
		 * Annotate each branch by updating its data field
		 * according to scores and user's annotation preferences. 
		 */
		NodeData nd = null;
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;
			
			if (!nodeDataList.containsKey(node)) {
				node.setData(null);
				continue;
			}
			nd = nodeDataList.get(node);
				
			double f1 = nd.mainfreq;
			double f2 = nd.alt1freqs;
			double f3 = nd.alt2freqs;
			double quarc = nd.quartcount;
			double effni = nd.effn + 0.0;
			
			if ( Math.abs((f1+f2+f3) - effni) > 0.01 ) {
				//System.err.println("Adjusting effective N from\t" + effni + "\tto\t" + (f1 + f2 + f3) + ". This should only happen as a result of polytomies in gene trees.");
				effni = f1 + f2 + f3;
			}
			
			//Long sum = p+a1+a2;
			
			Posterior post = new Posterior(
					f1,f2,f3,(double)effni, options.getLambda());
			Posterior post2 = new Posterior(
					f2,f1,f3, (double)effni, options.getLambda());
			Posterior post3 = new Posterior(
					f3,f1,f2, (double)effni, options.getLambda());
			
			nd.postQ1 = post;
			nd.postQ2 = post2;
			nd.postQ3 = post3;
			
			double bl = post.branchLength();
			
			node.setParentDistance(bl);
			nd.setString(this.getBranchAnnotation());		
			System.err.println(nd.data);
			node.setData(nd.toString2(this.getBranchAnnotation()));
				//i++;
		}
		
		return ret;
	}


	@Override
	Long getTotalCost(Vertex all) {
		System.err.println("Normalized score (portion of input quartet trees satisfied): " + 
				all._max_score/4./this.maxpossible);
		return (long) (all._max_score/4l);
	}


	@Override
	AbstractComputeMinCostTask<Tripartition> newComputeMinCostTask(AbstractInference<Tripartition> dlInference,
			Vertex all, IClusterCollection clusters) {
		return new WQComputeMinCostTask( (WQInference) dlInference, all,  clusters);
	}

	IClusterCollection newClusterCollection() {
		return new WQClusterCollection(GlobalMaps.taxonIdentifier.taxonCount());
	}
	
	WQDataCollection newCounter(IClusterCollection clusters) {
		return new WQDataCollection((WQClusterCollection)clusters, this);
	}



	@Override
	AbstractWeightCalculator<Tripartition> newWeightCalculator() {
		return new WQWeightCalculator(this);
	}


	@Override
	void setupMisc() {
		this.maxpossible = this.calculateMaxPossible();
		System.err.println("Number of quartet trees in the gene trees: " +
				this.maxpossible);
		
	}

	/**
	 * obsolete
	 */
	private void automaticallyDecideAlgorithm(int geneTreeTripartitonCountSize, int k){
		if (this.forceAlg != -1) {
			return;
		}
		if (k <= 0 || geneTreeTripartitonCountSize <= 0) {
			throw new RuntimeException("gene tree tripartition size or k not set properly");
		}
		if (this.forceAlg == -1) {
			this.forceAlg = ( GlobalMaps.taxonIdentifier.taxonCount() <= 32 || (geneTreeTripartitonCountSize < k*6)) ? 2 : 1;
		} else {
			throw new RuntimeException("Algorithm already set");
		}
	}
}
