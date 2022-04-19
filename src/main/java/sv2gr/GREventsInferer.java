package sv2gr;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import basic.Pair;
import population.popHierarchy.PopulationStructureFileReader;
import population.sv.utils.SimpleSVLocus;
import population.sv.utils.SimpleSVType;
import sv2gr.gr.GREvent;
import sv2gr.tree.RegionalTree;

/**
 * 
 * @author tanxu
 *
 */
public class GREventsInferer {
	/**
	 * regional trees with assigned SVs
	 */
	private final List<RegionalTree> regionalTrees;
	
	/**
	 * 
	 */
	private final PopulationStructureFileReader populationStructureFileReader;
	
	/**
	 * if true, treat all sample with HETER type genotype as Presence of the SV
	 * if false, all such samples will be treated as missing data
	 */
	private final boolean toTreatHeterGenotypeAsPresence;
	
	/**
	 * if true, the outgroup sample's genotype is always considered in the inference of GR event; 
	 * if false, the outgroup sample's genotype is only used as tie-breaker if there is exactly one single ingroup sample in sub-tree for both Presence and Absence genotype
	 * 
	 */
	private final boolean alwaysConsiderOutgroupSampleGenotype;
	
//	/**
//	 * a value between [0,1] indicate the minimal value of branch bootstrap value
//	 */
//	private final double minBootstrapValue;
	
	
//	/**
//	 * 
//	 */
//	private final boolean allowingPartiallyOverlappingGREventPath;
//	
//	/**
//	 * if true, check an additional condition that 
//	 * 		assume the GR event pairs occurred at time point t1 and t2,
//	 * 		only if there is no other GR event within the time period between t1 and t2 , the GR event pair will be considered to be retrieved
//	 * 
//	 * if false, simply retrieve all GR event pairs that meet the other conditions
//	 */
//	private final boolean onlyExtractContinuousGREventPairs;
	
	
//	private final String outgroupSampleName;
	
	//////////////////////////////
	/**
	 * 
	 */
	private Map<RegionalTree, List<Pair<GREvent, GREvent>>> treeTemporallyOrderedGREventPairsMap;
	
	
	
	public GREventsInferer(
			List<RegionalTree> regionalTrees, 
			PopulationStructureFileReader populationStructureFileReader,
//			double minBootstrapValue,
//			boolean allowingPartiallyOverlappingGREventPath,
			boolean toTreatHeterGenotypeAsPresence,
			boolean alwaysConsiderOutgroupSampleGenotype
//			boolean onlyExtractContinuousGREventPairs
//			, String outgroupSampleName
			) {
		super();
		this.regionalTrees = regionalTrees;
		this.populationStructureFileReader=populationStructureFileReader;
//		this.minBootstrapValue = minBootstrapValue;
//		this.allowingPartiallyOverlappingGREventPath = allowingPartiallyOverlappingGREventPath;
		this.toTreatHeterGenotypeAsPresence=toTreatHeterGenotypeAsPresence;
		this.alwaysConsiderOutgroupSampleGenotype=alwaysConsiderOutgroupSampleGenotype;
//		this.onlyExtractContinuousGREventPairs=onlyExtractContinuousGREventPairs;
//		this.outgroupSampleName = outgroupSampleName;
		
		//////////////////////////////
		this.inferGREvents();
//		this.resolveTemporalOrder();
	}
	
	void inferGREvents() {
		for(RegionalTree tree:this.regionalTrees) {
			tree.runPipelineToInferAllGREvents(this.populationStructureFileReader, this.toTreatHeterGenotypeAsPresence, this.alwaysConsiderOutgroupSampleGenotype);
		}
	}
	
//	/**
//	 * for every pair of GREvents g1 and g2, if the following conditions are all met
//	 * 
//	 * 1. g1 and g2 are on the same leaf-root path
//	 * 
//	 * 2. if g1 and g2 are overlapping, g1 and g2 must not have the same set of assigned branches
//	 * 
//	 * 3. g1 and g2 does not share any branch if {@link #allowingPartiallyOverlappingGREventPath} is false;
//	 * 
//	 * 4. tree path of g1 and g2 are not exact the same as the path from root node of ingroup samples to the root node of the tree!!!!
//	 * 		otherwise, the GR events occured before the divergence of the ingroup samples!!!!
//	 * 
//	 * 5. tree path of g1 and g2 are not exact the same as the path from outgroup leaf node to root node of the tree
//	 * 		otherwise, the GR events occured in the outgroup after its divergence with the ingroup samples!!!!
//	 * 
//	 * 6. tree path of g1 and g2 not all have bootstrap value < {@link #minBootstrapValue}
//	 * 		or at least one branch of g1 and one branch of g2 have bootstrap value >= {@link #minBootstrapValue}
//	 * 		note that for leaf node, its branch to parent node is always assigned bootstrap value 1!!!!
//	 * 		
//	 * then the GREvent with assigned path closer to the root node is considered occurred before the other one
//	 * 
//	 */
//	void resolveTemporalOrder() {
//		this.treeTemporallyOrderedGREventPairsMap=new HashMap<>();
//		
//		for(RegionalTree tree:this.regionalTrees) {
//			for(SimpleSVLocus sv1:tree.getSvLocusGREventMap().keySet()) {
//				
//				GREvent gr1=tree.getSvLocusGREventMap().get(sv1);
//				
//				for(SimpleSVLocus sv2:tree.getSvLocusGREventMap().keySet()) {
//					GREvent gr2=tree.getSvLocusGREventMap().get(sv2);
//					
//					if(gr1!=gr2) {
//						if(!DescendantToAncestralTreePathUtils.onSameLeafRootPath(gr1.getAssignedTreePath(), gr2.getAssignedTreePath())) {
//							//violate condition 1, skip
//							continue;
//						}
//						
//						if(DescendantToAncestralTreePathUtils.sameTreePath(gr1.getAssignedTreePath(), gr2.getAssignedTreePath())) {
//							//violate condition 2, skip
//							continue;
//						}
//						
//						if(!this.allowingPartiallyOverlappingGREventPath 
//								&& DescendantToAncestralTreePathUtils.getSharedTreePath(gr1.getAssignedTreePath(), gr2.getAssignedTreePath()).getBranchNumOnPath()>0) {
//							//violate condition 3, skip
//							continue;
//						}
//						//////////
//						if(DescendantToAncestralTreePathUtils.sameTreePath(gr1.getAssignedTreePath(), tree.getTree().getPathFromIngroupRootNodeToRootNode())
//								||DescendantToAncestralTreePathUtils.sameTreePath(gr2.getAssignedTreePath(), tree.getTree().getPathFromIngroupRootNodeToRootNode())) {
//							//violate condition 4, skip
//							continue;
//						}
//						
//						if(DescendantToAncestralTreePathUtils.sameTreePath(gr1.getAssignedTreePath(), tree.getTree().getPathFromOutgroupLeafNodeToRootNode())
//								||DescendantToAncestralTreePathUtils.sameTreePath(gr2.getAssignedTreePath(), tree.getTree().getPathFromOutgroupLeafNodeToRootNode())) {
//							//violate condition 5, skip
//							continue;
//						}
//						
//						
//						///////////check bootstrap value
//						boolean gr1TreePathBranchAllHaveValidBoostrap=false;
//						boolean gr2TreePathBranchAllHaveValidBoostrap=false;
//						for(int i=0;i<gr1.getAssignedTreePath().getNodesFromDescendantToAncestor().size()-1;i++) {
//							if(gr1.getAssignedTreePath().getNodesFromDescendantToAncestor().get(i).getBootstrap()>=this.minBootstrapValue) {
//								gr1TreePathBranchAllHaveValidBoostrap=true;
//							}
//						}
//						for(int i=0;i<gr2.getAssignedTreePath().getNodesFromDescendantToAncestor().size()-1;i++) {
//							TreeNode node=gr2.getAssignedTreePath().getNodesFromDescendantToAncestor().get(i);
//							Double bs=node.getBootstrap();
//							if(bs==null) {
//								System.out.println();
//							}
//							if(gr2.getAssignedTreePath().getNodesFromDescendantToAncestor().get(i).getBootstrap()>=this.minBootstrapValue) {
//								gr2TreePathBranchAllHaveValidBoostrap=true;
//							}
//						}
//						if(!(gr1TreePathBranchAllHaveValidBoostrap&&gr2TreePathBranchAllHaveValidBoostrap)) {
//							//violate condition 6, skip
//							continue;
//						}
//						
//						//////////////////
//						if(this.onlyExtractContinuousGREventPairs) {
//							DescendantToAncestralTreePath fullPath=DescendantToAncestralTreePathUtils.getShortestCoveringTreePath(gr1.getAssignedTreePath(), gr2.getAssignedTreePath());
//							
//							boolean otherGREventPathBetweenGR1AndGR2Found=false;
//							for(GREvent gr:tree.getSvLocusGREventMap().values()) {
//								if(gr!=gr1 && gr!=gr2) {
//									if(DescendantToAncestralTreePathUtils.p1IsFullyInsideOFP2(gr.getAssignedTreePath(), fullPath)) {
//										otherGREventPathBetweenGR1AndGR2Found=true;
//									}
//								}
//							}
//							
//							if(otherGREventPathBetweenGR1AndGR2Found)
//								continue;
//						}
//						
//						
//						
//						/////////////////////////
//						if(gr1.getAssignedTreePath().getAncestralNode().getAllAncestralNodes().contains(gr2.getAssignedTreePath().getAncestralNode())) {
//							//gr1's tree path is at the downstream (descendant) of gr2
//							//gr2 occurred before gr1
//							Pair<GREvent, GREvent> pair=new Pair<>(gr2, gr1);
//							
//							if(!this.treeTemporallyOrderedGREventPairsMap.containsKey(tree))
//								this.treeTemporallyOrderedGREventPairsMap.put(tree, new ArrayList<>());
//							
//							if(!this.treeTemporallyOrderedGREventPairsMap.get(tree).contains(pair))
//								this.treeTemporallyOrderedGREventPairsMap.get(tree).add(pair);
//							
//						}else {
//							//gr2's tree path is at the downstream (descendant) of gr1
//							//gr1 occurred before gr2
//							Pair<GREvent, GREvent> pair=new Pair<>(gr1, gr2);
//							
//							if(!this.treeTemporallyOrderedGREventPairsMap.containsKey(tree))
//								this.treeTemporallyOrderedGREventPairsMap.put(tree, new ArrayList<>());
//							
//							if(!this.treeTemporallyOrderedGREventPairsMap.get(tree).contains(pair))
//								this.treeTemporallyOrderedGREventPairsMap.get(tree).add(pair);
//						}
//					}
//				}
//			}
//		}
//	}

	/**
	 * @return the treeTemporallyOrderedGREventPairsMap
	 */
	public Map<RegionalTree, List<Pair<GREvent, GREvent>>> getTreeTemporallyOrderedGREventPairsMap() {
		return treeTemporallyOrderedGREventPairsMap;
	}

	/**
	 * @return the regionalTrees
	 */
	public List<RegionalTree> getRegionalTrees() {
		return regionalTrees;
	}
	
	
	/**
	 * 
	 * @return 
	 */
	public Map<SimpleSVType, Pair<Integer,Integer>> getSVTypeNumWithInferedGREventMap(){
		Map<SimpleSVType, Integer> svTypeReversedGREventNumMap = new HashMap<>();
		Map<SimpleSVType, Integer> svTypeNonReversedGREventNumMap = new HashMap<>();
		Set<SimpleSVType> includedTypes=new HashSet<>();
		for(RegionalTree tree:this.regionalTrees) {
			for(SimpleSVLocus sv:tree.getSvLocusGREventMap().keySet()) {
				GREvent gr=tree.getSvLocusGREventMap().get(sv);
				includedTypes.add(sv.getType());
				if(gr.isReversalOfSV()) {
					if(!svTypeReversedGREventNumMap.containsKey(sv.getType())) {
						svTypeReversedGREventNumMap.put(sv.getType(), 0);
					}
					svTypeReversedGREventNumMap.put(sv.getType(), svTypeReversedGREventNumMap.get(sv.getType())+1);
				}else {
					if(!svTypeNonReversedGREventNumMap.containsKey(sv.getType())) {
						svTypeNonReversedGREventNumMap.put(sv.getType(), 0);
					}
					svTypeNonReversedGREventNumMap.put(sv.getType(), svTypeNonReversedGREventNumMap.get(sv.getType())+1);
				}
			}
		}
		
		
		//////////////
		Map<SimpleSVType, Pair<Integer,Integer>> ret =new HashMap<>();
		
		for(SimpleSVType type:includedTypes) {
			int reversed=0;
			int nonReversed=0;
			if(svTypeReversedGREventNumMap.containsKey(type)) {
				reversed=svTypeReversedGREventNumMap.get(type);
			}
			if(svTypeNonReversedGREventNumMap.containsKey(type)) {
				nonReversed=svTypeNonReversedGREventNumMap.get(type);
			}
			
			ret.put(type, new Pair<>(reversed, nonReversed));
		}
		
		return ret;
	}
	
	
}
