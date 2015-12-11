package phylonet.coalescent;
import phylonet.coalescent.BipartitionWeightCalculator.Results;
import phylonet.coalescent.Posterior;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Stack;

import phylonet.coalescent.BipartitionWeightCalculator.Quadrapartition;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;

public class WQInference extends AbstractInference<Tripartition> {
	
	int forceAlg = -1;
	long maxpossible;
	public WQInference(boolean rooted, boolean extrarooted, List<Tree> trees,
			List<Tree> extraTrees, boolean exactSolution, boolean duploss, int alg, int addExtra,
			boolean outputCompletedGenes, boolean outSearch, boolean run) {
		super(rooted, extrarooted, trees, extraTrees, exactSolution, 
				addExtra, outputCompletedGenes, outSearch, run);
		
		this.forceAlg = alg;
	}

	
	public void scoreGeneTree(Tree st) {

		mapNames();

		IClusterCollection clusters = newClusterCollection();


		this.dataCollection = newCounter(clusters);
		weightCalculator = newWeightCalculator();

		WQDataCollection wqDataCollection = (WQDataCollection) this.dataCollection;
		wqDataCollection.preProcess(this);
		wqDataCollection.initializeWeightCalculator(this);
		this.maxpossible = wqDataCollection.calculateMaxPossible();
		//System.err.println("Number of quartet trees in the gene trees: "+this.maxpossible);

		System.err.println(this.maxpossible);
		
		Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
		long sum = 0l;
		long maxsum = 0l;
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
				for (int i = 0; i < childbslist.size(); i++) {
					for (int j = i+1; j < childbslist.size(); j++) {
						for (int k = j+1; k < childbslist.size(); k++) {
							Tripartition trip = new Tripartition(childbslist.get(i),  childbslist.get(j), childbslist.get(k));
							Long s = weightCalculator.getWeight(trip, null);
							sum += s;
							//Long m = this.dataCollection.maxPossibleScore(trip);
							//((STINode)node).setData(s*100l/m);
							//maxsum += m;
						}
					}					       
				}
				if (childbslist.size() > 3) {
					for (STITreeCluster chid :childbslist) {
						System.err.print(chid.getClusterSize()+" ");
					}
					System.err.println(" (polytomy)");
				}
			}
		}
/*		if (4l*this.maxpossible != maxsum) {
			throw new RuntimeException("Hmm... "+maxsum+" "+4l*this.maxpossible);
		}*/
		System.err.println("Quartet score is: " + sum/4l);
		System.err.println("Normalized quartet score is: "+ (sum/4l+0.)/this.maxpossible);
		//System.out.println(st.toNewickWD());
		this.scoreBranches(st);
	}
	
	
	
	public void scoreBranches(Tree st) {

		weightCalculator = new BipartitionWeightCalculator(this);
		
		BipartitionWeightCalculator weightCalculator2 = (BipartitionWeightCalculator) weightCalculator;
		WQDataCollection wqDataCollection = (WQDataCollection) this.dataCollection;
		//wqDataCollection.initializeWeightCalculator(this);
		
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
		List<Long> mainfreqs = new ArrayList<Long>();
		List<Long> alt1freqs = new ArrayList<Long>();
		List<Long> alt2freqs = new ArrayList<Long>();
		List<Long> quartcount = new ArrayList<Long>();
		List<Integer> effn = new ArrayList<Integer>();
		
		
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;
			if (node.isLeaf()) {
				stack.push((STITreeCluster) node.getData());
			} else {
				STITreeCluster cluster = (STITreeCluster) node.getData();				
				if (node.isRoot() || node.getChildCount() > 2 || (node.getParent().getChildCount() >3) ||
						(node.getParent().getChildCount() >2 && !node.getParent().isRoot()) ) {
					for (int i =0; i< node.getChildCount(); i++) {
						stack.pop();
					}
					stack.push(cluster);
					continue;
				}
				
				STITreeCluster c1 = stack.pop();
				STITreeCluster c2 = stack.pop();
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
				} else if (node.getParent().isRoot() && node.getParent().getChildCount() == 2) {
					continue;
				}
				else {
					remaining = ((STITreeCluster)node.getParent().getData()).complementaryCluster();
				}
				Quadrapartition quadm = weightCalculator2.new Quadrapartition
						(c1,  c2, sister, remaining);
				STBipartition bmain = new STBipartition(cluster, cluster.complementaryCluster());
				
				Results s = weightCalculator2.getWeight(quadm);
				Long p = s.qs;
				mainfreqs.add(s.qs);
				effn.add(s.effn);
				
				
				Quadrapartition quad2 = weightCalculator2.new Quadrapartition
						(c1, sister, c2, remaining);
				STITreeCluster c1plussis = new STITreeCluster();
				c1plussis.setCluster((BitSet) c1.getBitSet().clone());
				c1plussis.getBitSet().or(sister.getBitSet());
				STBipartition b2 = new STBipartition(c1plussis, c1plussis.complementaryCluster());
				s = weightCalculator2.getWeight(quad2);
				Long a1 = s.qs;

				alt1freqs.add(a1);
				
				Quadrapartition quad3 = weightCalculator2.new Quadrapartition
						(c1, remaining, c2, sister);
				STITreeCluster c1plusrem = new STITreeCluster();
				c1plusrem.setCluster((BitSet) c1.getBitSet().clone());
				c1plusrem.getBitSet().or(remaining.getBitSet());
				STBipartition b3 = new STBipartition(c1plusrem, c1plusrem.complementaryCluster());
				s = weightCalculator2.getWeight(quad3);
				Long a2 = s.qs;
				alt2freqs.add(a2);


				
				quartcount.add( (c1.getClusterSize()+0l)
						* (c2.getClusterSize()+0l)
						* (sister.getClusterSize()+0l)
						* (remaining.getClusterSize()+0l));
				

				Posterior pst_tmp = new Posterior((double)p,(double)a1,(double)a2,(double)s.effn);
				double post_m = pst_tmp.getPost();
				pst_tmp = new Posterior((double)a1,(double)p,(double)a2,(double)s.effn);

				double post_a1 = pst_tmp.getPost();
				//pst_tmp =  new Posterior((double)a2,(double)p,(double)a1,(double)numTrees);
				double post_a2 = Math.max(0.,1.0 - post_a1 - post_m);
				
				
				System.err.println(quadm +
						" [" + bmain.toString2() +"] : "+post_m +" ** f1 = "+p+
						" f2 = "+a1+" f3 = "+a2+" effective_n = "+ (double)s.effn+" **");
			
				if (this.getBranchAnnotation() == 4){
					
						System.err.println(quad2 +
								" ["+b2.toString2()+"] : "+post_a1+ " ** f1 = "+a1+
								" f2 = "+p+" f3 = "+a2+" effective_n = "+ (double)s.effn+" **");
						System.err.println(quad3 +
								" ["+b3.toString2()+"] : "+post_a2+ " ** f1 = "+a2+
								" f2 = "+p+" f3 = "+a1+" effective_n = "+ (double)s.effn+" **");
					
				}
				//System.err.println(quad2+" : "+post_a1);
				//System.err.println(quad3+" : "+post_a2);
			}
		}
		int i = 0;
		for (TNode n: st.postTraverse()) {
			STINode node = (STINode) n;
			if (node.isLeaf() || node.isRoot() ||
					(node.getParent().isRoot() && node.getParent().getChildCount() ==2) ||
					  node.getChildCount() > 2 || (node.getParent().getChildCount() >3) ||
						(node.getParent().getChildCount() >2 && !node.getParent().isRoot())) {
				node.setData(null);
			} else{
				Long p = mainfreqs.get(i);
				Long a1 = alt1freqs.get(i);
				Long a2 = alt2freqs.get(i);
				Long quarc = quartcount.get(i);
				Integer effni = effn.get(i);
				Long sum = p+a1+a2;
				
				Double bl = -Math.log(1.5*(1.0-((p+.0)/sum)));
				if (bl.isInfinite()) {
					bl = 10.;
				}
				node.setParentDistance(bl);
				if (this.getBranchAnnotation() == 0){
					node.setData(null);
				} else if (this.getBranchAnnotation() == 1){
					node.setData((p+.0)/sum*100);
				} else {
					Posterior pst_tmp = new Posterior((double)p,(double)a1,(double)a2,(double)effni);
					double post_m = pst_tmp.getPost();
					 
					if (this.getBranchAnnotation() == 3) {
						node.setData(post_m);
					} else if (this.getBranchAnnotation() == 2) {
						pst_tmp = new Posterior((double)a1,(double)p,(double)a2,(double)effni);
						double post_a1 = pst_tmp.getPost();
						//pst_tmp =  new Posterior((double)a2,(double)p,(double)a1,(double)numTrees);
						double post_a2 = 1.0 - post_a1 - post_m;
						
						node.setData("[q1="+(p+.0)/sum+";q2="+(a1+.0)/sum+";q3="+(a2+.0)/sum+
									 ";f1="+p+";f2="+a1+";f3="+a2+
									 ";pp1="+post_m+";pp2="+post_a1+";pp3="+post_a2+
									 ";QC="+quarc+"]");
					}  else if (this.getBranchAnnotation() == 4) {
						pst_tmp = new Posterior((double)a1,(double)p,(double)a2,(double)effni);
						double post_a1 = pst_tmp.getPost();
						//pst_tmp =  new Posterior((double)a2,(double)p,(double)a1,(double)numTrees);
						double post_a2 = 1.0 - post_a1 - post_m;
						
						node.setData("[q1="+(p+.0)/sum+";q2="+(a1+.0)/sum+";q3="+(a2+.0)/sum+
									 ";f1="+p+";f2="+a1+";f3="+a2+
									 ";pp1="+post_m+";pp2="+post_a1+";pp3="+post_a2+
									 ";QC="+quarc+"]");
					} 
				}
				i++;
			} 
		}
		
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
		return new WQDataCollection((WQClusterCollection)clusters, this.forceAlg, this);
	}



	@Override
	AbstractWeightCalculator<Tripartition> newWeightCalculator() {
		return new WQWeightCalculator(this);
	}

}
