package sv2gr.tree;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.google.common.base.Objects;

import population.popHierarchy.PopulationStructureFileReader;
import population.sv.utils.SimpleSVLocus;
import population.utils.Genotype;
import sv2gr.SvToRegionalTreeAssigner;
import sv2gr.gr.GREvent;

/**
 * 
 * @author tanxu
 *
 */
public class RegionalTree {
	/**
	 * regional tree ID;
	 */
	private final int ID;
	/**
	 * 
	 */
	private final int windowSize;
	/**
	 * 
	 */
	private final String chrom;
	/**
	 * 
	 */
	private final int start;
	/**
	 * 
	 */
	private final int end;
	
	/**
	 * 
	 */
	private final Tree tree;
	
	/////////////////////////////////
	/**
	 * all {@link SimpleSVLocus}s assigned to this {@link RegionalTree}
	 * see {@link SvToRegionalTreeAssigner} for more details
	 */
	private List<SimpleSVLocus> coveredSVLocusList;
	
//	private Map<String, TreeNode> sampleNameTreeLeafNodeMapMap;
//	private String outgroupSampleName;
	
	/**
	 * map from {@link SimpleSVLocus} to the inferred {@link GREvent} of the SV
	 * note that only {@link GREvent} assigned to tree path with at least one branch within the clade of the ingroup samples are stored 
	 * 		for {@link GREvent} inferred to be located on 
	 * 			1. tree branch from outgroup sample leaf node to root node or
	 * 			2. tree branch from root node of clade of all ingroup sample to root node of the tree
	 * 		they will not be stored!!!!!
	 * 		
	 * map value can be null for some {@link SimpleSVLocus}s
	 */
	private Map<SimpleSVLocus, GREvent> svLocusGREventMap;
	
	/**
	 * 
	 * @param ID
	 * @param windowSize
	 * @param chrom
	 * @param start
	 * @param end
	 * @param tree
	 */
	public RegionalTree(int ID, int windowSize, String chrom, int start, int end, Tree tree) {
		super();
		this.ID=ID;
		this.windowSize=windowSize;
		this.chrom = chrom;
		this.start = start;
		this.end = end;
		this.tree = tree;
		
		
		/////////////////////
		this.coveredSVLocusList=new ArrayList<>();
	}
	
	/**
	 * add SimpleSVLocus that is determined to be covered by this regional tree
	 * @param locus
	 */
	public void addCoveredSVLocus(SimpleSVLocus locus) {
		this.coveredSVLocusList.add(locus);
	}
	
	
	///////////////////////////////////
	/**
	 * run the full pipline to TODO
	 * 
	 * must be invoked after all {@link SimpleSVLocus} covered by this regional tree are added!
	 * @param reader
	 * @param toTreatHeterGenotypeAsPresence if true, treat all sample with HETER type genotype as Presence of the SV; if false, all such samples will be treated as missing data
	 * @param alwaysConsiderOutgroupSampleGenotype if true, the outgroup sample genotype is always included in the inference of GR event; if false, it is only used as tie-breaker if there is only one single ingroup sample in sub-tree for both Presence and Absence genotype
	 */
	public void runPipelineToInferAllGREvents(PopulationStructureFileReader reader, boolean toTreatHeterGenotypeAsPresence, boolean alwaysConsiderOutgroupSampleGenotype) {
		///
		this.recodeSampleGenotypes(toTreatHeterGenotypeAsPresence);
		
		////
		this.assignSampleIndexToLeafNodes(reader.getSampleNameIndexMap());
		
		this.setSVLocusForAllNodes();
		
		this.assignSVGenotypesToAllNodes();
		
		this.identifyTreeNodesInSubtreeForEachSV();
		
//		this.inferGREventForEachSV(reader, alwaysConsiderOutgroupSampleGenotype);
		this.inferGREventForEachSV2(reader, alwaysConsiderOutgroupSampleGenotype);
	}
	
	
	/////////////////////////////////preprocessing step : recode ingroup sample genotype as required
	//
	void recodeSampleGenotypes(boolean toTreatHeterGenotypeAsPresence) {
		for(SimpleSVLocus sv:this.coveredSVLocusList) {
			for(int sampleIndex:sv.getSampleIndexGenotypeMap().keySet()) {
				Genotype gt=sv.getSampleIndexGenotypeMap().get(sampleIndex);
				
				if(gt.equals(Genotype.HETER)) {
					if(toTreatHeterGenotypeAsPresence) {
						sv.recodeSampleGenotype(sampleIndex, Genotype.PRESENCE);
					}else {
						sv.recodeSampleGenotype(sampleIndex, Genotype.MISSING);
					}
				}
			}
		}
	}
	
	/////////////////////////////////////step 1
	
	/**
	 * assign the sample index to each leaf node containing the corresponding sample;
	 * 
	 * the sample index should be consistent with the {@link SimpleSVLocus#getSampleIndexGenotypeMap()}
	 * @param sampleNameIndexMap the map contains a super set of all samples on tree
	 */
	void assignSampleIndexToLeafNodes(Map<String, Integer> sampleNameIndexMap) {
		for(String label: this.tree.getLeafLabelTreeNodeMap().keySet()) {
			if(!sampleNameIndexMap.containsKey(label)) {
				throw new IllegalArgumentException("leaf node sample is not found in given map!");
			}
			TreeNode leaf=this.tree.getLeafLabelTreeNodeMap().get(label);
			int index=sampleNameIndexMap.get(label);
			leaf.setSampleIndex(index);
		}
	}
	
	/**
	 */
	void setSVLocusForAllNodes() {
		//initialize the sv - genotype map for every tree node (including leaf and internal ones)
		for(TreeNode treeNode:this.tree.getNodeIDMap().values()) {
			treeNode.setCoveredSVLocusList(coveredSVLocusList);
		}
	}
	
	/**
	 * assign genotype of each SV to each tree nodes including both leaf nodes and internal nodes
	 */
	void assignSVGenotypesToAllNodes() {
		this.tree.getRootNode().inferSVGenotypes();
	}
	
	/**
	 * identify the tree nodes that should be included in the sub-tree for each SV for GR event inference
	 */
	void identifyTreeNodesInSubtreeForEachSV() {
		this.tree.getRootNode().inferInSubTree();
	}
	
	
	/**
	 * for each SV covered by this Tree, infer the GR event if possible
	 * //TODO testing
	 * @param alwaysConsiderOutgroupSampleGenotype if true, the outgroup sample genotype is always included in the inference of GR event; if false, it is only used as tie-breaker if there is only one single ingroup sample in sub-tree for both Presence and Absence genotype
	 */
	void inferGREventForEachSV(PopulationStructureFileReader reader, boolean alwaysConsiderOutgroupSampleGenotype) {
		this.svLocusGREventMap=new LinkedHashMap<>();
		
		for(SimpleSVLocus sv:this.coveredSVLocusList) {
			//find out all leaf nodes in subtree for each possible genotype; these leaf nodes on subtree may be non-leaf node on the original tree!!!!
			//notes that only Presence and absence genotypes can be found for leaf nodes in subtree
			Map<Genotype, List<TreeNode>> genotypeLeafNodesOnSubTreeMap = new HashMap<>();
			
			for(int id:this.tree.getNodeIDMap().keySet()) {
				TreeNode node = this.tree.getNodeIDMap().get(id);
				
				if(node.getSvLocusNodeInSubTreeMap().get(sv)) {//included in subtree
					if(node.isLeaf()) {//leaf node on original tree, always be included, trivial
						Genotype gt=node.getSvLocusGenotypeMap().get(sv);
						if(!genotypeLeafNodesOnSubTreeMap.containsKey(gt)) {
							genotypeLeafNodesOnSubTreeMap.put(gt, new ArrayList<>());
						}
						genotypeLeafNodesOnSubTreeMap.get(gt).add(node);
					}else {//non-leaf node on the original tree, check if it is leaf node on the subtree 
						//
						boolean childNodeIncluded=false;
						for(TreeNode child:node.getChildNodeList()) {
							if(child.getSvLocusNodeInSubTreeMap().get(sv)) {
								childNodeIncluded=true;
							}
						}
						if(!childNodeIncluded) {
							Genotype gt=node.getSvLocusGenotypeMap().get(sv);
							if(!genotypeLeafNodesOnSubTreeMap.containsKey(gt)) {
								genotypeLeafNodesOnSubTreeMap.put(gt, new ArrayList<>());
							}
							genotypeLeafNodesOnSubTreeMap.get(gt).add(node);
						}
					}
				}
			}
			
			
			if(alwaysConsiderOutgroupSampleGenotype) {
				/////////////////////for the two genotypes (0/0 and 1/1), 
				///////////////////////check if one genotype has 1 leaf node on subtree and the other genotype has 2 or more leaf node on tree
				int absenceGenotypeNum=genotypeLeafNodesOnSubTreeMap.containsKey(Genotype.ABSENCE)?genotypeLeafNodesOnSubTreeMap.get(Genotype.ABSENCE).size():0;
				int presenceGenotypeNum=genotypeLeafNodesOnSubTreeMap.containsKey(Genotype.PRESENCE)?genotypeLeafNodesOnSubTreeMap.get(Genotype.PRESENCE).size():0;
				
				if(absenceGenotypeNum==1 && presenceGenotypeNum>=2) {//GR event of the reverse of the SV occurred on the tree path starting from the node with absence genotype (0/0)
					//identify the tree branch path
					TreeNode start = genotypeLeafNodesOnSubTreeMap.get(Genotype.ABSENCE).get(0);
					List<TreeNode> treeNodesOnPath = new ArrayList<>();
					treeNodesOnPath.add(start);
					
					//trace back the path until the first ancestral node with one leaf node on sub-tree with genotype 'presence' (1/1) is encountered or the root node is reached
					TreeNode ancestralNode=start.getParent();
					boolean ancestralNodeHasLeafNodeOnSubtreeWithPresenceGenotype=false;
					if(ancestralNode!=null) {
						treeNodesOnPath.add(ancestralNode);
						for(TreeNode child:ancestralNode.getChildNodeList()) {
							if(child.getSvLocusNodeInSubTreeMap().get(sv) && Objects.equal(child.getSvLocusGenotypeMap().get(sv), Genotype.PRESENCE)) {
								ancestralNodeHasLeafNodeOnSubtreeWithPresenceGenotype=true;
							}
						}
					}
					while(!ancestralNodeHasLeafNodeOnSubtreeWithPresenceGenotype && ancestralNode!=null) {//neither the ancestral node with the one leaf node on sub-tree with genotype 'presence' (1/1) is encountered nor the root node is reached
						ancestralNode=ancestralNode.getParent();
						ancestralNodeHasLeafNodeOnSubtreeWithPresenceGenotype=false;
						if(ancestralNode!=null) {
							treeNodesOnPath.add(ancestralNode);
							for(TreeNode child:ancestralNode.getChildNodeList()) {
								if(child.getSvLocusNodeInSubTreeMap().get(sv) && Objects.equal(child.getSvLocusGenotypeMap().get(sv), Genotype.PRESENCE)) {
									ancestralNodeHasLeafNodeOnSubtreeWithPresenceGenotype=true;
								}
							}
						}
					}
					
					//////////////////////////////////
					DescendantToAncestralTreePath treePath=new DescendantToAncestralTreePath(treeNodesOnPath);
					
					if(treePath.getBranchNumOnPath()>0) {//only add the GREvent if at least one branch is included
						GREvent grEvent=new GREvent(sv, this, true, treePath);
						
						if(DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromIngroupRootNodeToRootNode())
								||DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromOutgroupLeafNodeToRootNode())) {
							//skip if inferred GR event is only located on branches outside of the ingroup clade
						}else {
							this.svLocusGREventMap.put(sv, grEvent);
						}
					}
	
					
				}else if(absenceGenotypeNum>=2 && presenceGenotypeNum==1) {//GR event consistent with the SV occurred on the tree path starting from the node with presence genotype (0/0)
					//identify the tree branch path
					TreeNode start = genotypeLeafNodesOnSubTreeMap.get(Genotype.PRESENCE).get(0);
					List<TreeNode> treeNodesOnPath = new ArrayList<>();
					treeNodesOnPath.add(start);
					
					//trace back the path until the first ancestral node with one leaf node on sub-tree with genotype 'presence' (1/1) is encountered or the root node is reached
					TreeNode ancestralNode=start.getParent();
					boolean ancestralNodeHasLeafNodeOnSubtreeWithAbsenceGenotype=false;
					if(ancestralNode!=null) {
						treeNodesOnPath.add(ancestralNode);
						for(TreeNode child:ancestralNode.getChildNodeList()) {
							if(child.getSvLocusNodeInSubTreeMap().get(sv) && Objects.equal(child.getSvLocusGenotypeMap().get(sv), Genotype.ABSENCE)) {
								ancestralNodeHasLeafNodeOnSubtreeWithAbsenceGenotype=true;
							}
						}
					}
					while(!ancestralNodeHasLeafNodeOnSubtreeWithAbsenceGenotype && ancestralNode!=null) {//neither the ancestral node with the one leaf node on sub-tree with genotype 'presence' (1/1) is encountered nor the root node is reached
						ancestralNode=ancestralNode.getParent();
						ancestralNodeHasLeafNodeOnSubtreeWithAbsenceGenotype=false;
						if(ancestralNode!=null) {
							treeNodesOnPath.add(ancestralNode);
							for(TreeNode child:ancestralNode.getChildNodeList()) {
								if(child.getSvLocusNodeInSubTreeMap().get(sv) && Objects.equal(child.getSvLocusGenotypeMap().get(sv), Genotype.ABSENCE)) {
									ancestralNodeHasLeafNodeOnSubtreeWithAbsenceGenotype=true;
								}
							}
						}
					}
					//////////////////////////////////
					DescendantToAncestralTreePath treePath=new DescendantToAncestralTreePath(treeNodesOnPath);
					
					if(treePath.getBranchNumOnPath()>0) {//only add the GREvent if at least one branch is included
						GREvent grEvent=new GREvent(sv, this, false, treePath);
						
	
						if(DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromIngroupRootNodeToRootNode())
								||DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromOutgroupLeafNodeToRootNode())) {
							//skip if infered GR event is only located on branches outside of the ingroup clade
						}else {
							this.svLocusGREventMap.put(sv, grEvent);
						}
					}
				}else {//no GR event inferred
					//skip
				}
			}else {//only consider outgroup sample when necessary
				//find out the number of leaf node with with presence and absence genotype for all ingroup samples on the sub-tree
				//also find out the genotype of outgroup sample
				Genotype outgroupSampleGenotype=null;
				int ingroupSampleNumWithAbsenceGenotypeNum=0;
				int ingroupSampleNumWithPresenceGenotypeNum=0;
				
				if(genotypeLeafNodesOnSubTreeMap.containsKey(Genotype.ABSENCE)) {
					for(TreeNode node:genotypeLeafNodesOnSubTreeMap.get(Genotype.ABSENCE)) {
						//check if the tree node is within the clade containing all ingroup samples
						if(node.getParent()!=null) {//not root node
							//TODO
							if(reader.getAllOutgroupSampleIndices().contains(node.getSampleIndex())) {//outgroup sample node
								outgroupSampleGenotype=Genotype.ABSENCE;
							}else {//not root node, not outgroup sample node, must be within the clade of all ingroup samples
								ingroupSampleNumWithAbsenceGenotypeNum++;
							}
						}
					}
				}
				if(genotypeLeafNodesOnSubTreeMap.containsKey(Genotype.PRESENCE)) {
					for(TreeNode node:genotypeLeafNodesOnSubTreeMap.get(Genotype.PRESENCE)) {
						if(reader.getAllOutgroupSampleIndices().contains(node.getSampleIndex())) {//outgroup sample node
							outgroupSampleGenotype=Genotype.PRESENCE;
						}else {//not root node, not outgroup sample node, must be within the clade of all ingroup samples
							ingroupSampleNumWithPresenceGenotypeNum++;
						}
					}
				}
				/////////////////////////////////////
				if(ingroupSampleNumWithAbsenceGenotypeNum==1 && ingroupSampleNumWithPresenceGenotypeNum>=2 //
						|| (ingroupSampleNumWithAbsenceGenotypeNum==1 && ingroupSampleNumWithPresenceGenotypeNum==1)&&(outgroupSampleGenotype!=null && outgroupSampleGenotype.equals(Genotype.PRESENCE))
						////there is only one ingroup sample on the sub-tree for both presnce and absence genotype, check the outgroup sample genotype
						
						) {//GR event of the reverse of the SV occurred on the tree path starting from the node with absence genotype (0/0)
					//identify the tree branch path
					TreeNode start = genotypeLeafNodesOnSubTreeMap.get(Genotype.ABSENCE).get(0);
					List<TreeNode> treeNodesOnPath = new ArrayList<>();
					treeNodesOnPath.add(start);
					
					//trace back the path until the first ancestral node with one leaf node on sub-tree with genotype 'presence' (1/1) is encountered or the root node is reached
					TreeNode ancestralNode=start.getParent();
					boolean ancestralNodeHasLeafNodeOnSubtreeWithPresenceGenotype=false;
					if(ancestralNode!=null) {
						treeNodesOnPath.add(ancestralNode);
						for(TreeNode child:ancestralNode.getChildNodeList()) {
							if(child.getSvLocusNodeInSubTreeMap().get(sv) && Objects.equal(child.getSvLocusGenotypeMap().get(sv), Genotype.PRESENCE)) {
								ancestralNodeHasLeafNodeOnSubtreeWithPresenceGenotype=true;
							}
						}
					}
					while(!ancestralNodeHasLeafNodeOnSubtreeWithPresenceGenotype && ancestralNode!=null) {//neither the ancestral node with the one leaf node on sub-tree with genotype 'presence' (1/1) is encountered nor the root node is reached
						ancestralNode=ancestralNode.getParent();
						ancestralNodeHasLeafNodeOnSubtreeWithPresenceGenotype=false;
						if(ancestralNode!=null) {
							treeNodesOnPath.add(ancestralNode);
							for(TreeNode child:ancestralNode.getChildNodeList()) {
								if(child.getSvLocusNodeInSubTreeMap().get(sv) && Objects.equal(child.getSvLocusGenotypeMap().get(sv), Genotype.PRESENCE)) {
									ancestralNodeHasLeafNodeOnSubtreeWithPresenceGenotype=true;
								}
							}
						}
					}
					//////////////////////////////////
					DescendantToAncestralTreePath treePath=new DescendantToAncestralTreePath(treeNodesOnPath);
					
					if(treePath.getBranchNumOnPath()>0) {//only add the GREvent if at least one branch is included
						GREvent grEvent=new GREvent(sv, this, true, treePath);
						
						if(DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromIngroupRootNodeToRootNode())
								||DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromOutgroupLeafNodeToRootNode())) {
							//skip if infered GR event is only located on branches outside of the ingroup clade
						}else {
							this.svLocusGREventMap.put(sv, grEvent);
						}
					}
					
				}else if(ingroupSampleNumWithAbsenceGenotypeNum>=2 && ingroupSampleNumWithPresenceGenotypeNum==1
						|| (ingroupSampleNumWithAbsenceGenotypeNum==1 && ingroupSampleNumWithPresenceGenotypeNum==1)&&(outgroupSampleGenotype!=null && outgroupSampleGenotype.equals(Genotype.ABSENCE))
						
						) {//GR event consistent with the SV occurred on the tree path starting from the node with presence genotype (0/0)
					//identify the tree branch path
					TreeNode start = genotypeLeafNodesOnSubTreeMap.get(Genotype.PRESENCE).get(0);
					List<TreeNode> treeNodesOnPath = new ArrayList<>();
					treeNodesOnPath.add(start);
					
					//trace back the path until the first ancestral node with one leaf node on sub-tree with genotype 'presence' (1/1) is encountered or the root node is reached
					TreeNode ancestralNode=start.getParent();
					boolean ancestralNodeHasLeafNodeOnSubtreeWithAbsenceGenotype=false;
					if(ancestralNode!=null) {
						treeNodesOnPath.add(ancestralNode);
						for(TreeNode child:ancestralNode.getChildNodeList()) {
							if(child.getSvLocusNodeInSubTreeMap().get(sv) && Objects.equal(child.getSvLocusGenotypeMap().get(sv), Genotype.ABSENCE)) {
								ancestralNodeHasLeafNodeOnSubtreeWithAbsenceGenotype=true;
							}
						}
					}
					while(!ancestralNodeHasLeafNodeOnSubtreeWithAbsenceGenotype && ancestralNode!=null) {//neither the ancestral node with the one leaf node on sub-tree with genotype 'presence' (1/1) is encountered nor the root node is reached
						ancestralNode=ancestralNode.getParent();
						ancestralNodeHasLeafNodeOnSubtreeWithAbsenceGenotype=false;
						if(ancestralNode!=null) {
							treeNodesOnPath.add(ancestralNode);
							for(TreeNode child:ancestralNode.getChildNodeList()) {
								if(child.getSvLocusNodeInSubTreeMap().get(sv) && Objects.equal(child.getSvLocusGenotypeMap().get(sv), Genotype.ABSENCE)) {
									ancestralNodeHasLeafNodeOnSubtreeWithAbsenceGenotype=true;
								}
							}
						}
					}
					
					//////////////////////////////////
					DescendantToAncestralTreePath treePath=new DescendantToAncestralTreePath(treeNodesOnPath);
					
					if(treePath.getBranchNumOnPath()>0) {//only add the GREvent if at least one branch is included
						GREvent grEvent=new GREvent(sv, this, false, treePath);
						
	
						if(DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromIngroupRootNodeToRootNode())
								||DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromOutgroupLeafNodeToRootNode())) {
							//skip if infered GR event is only located on branches outside of the ingroup clade
						}else {
							this.svLocusGREventMap.put(sv, grEvent);
						}
					}
				}else {
					//skip
				}
			}
		}
	}
	
	/**
	 * for each SV covered by this Tree, infer the GR event if possible
	 * //TODO testing
	 * @param alwaysConsiderOutgroupSampleGenotype if true, the outgroup sample genotype is always included in the inference of GR event; if false, it is only used as tie-breaker if there is only one single ingroup sample in sub-tree for both Presence and Absence genotype
	 */
	void inferGREventForEachSV2(PopulationStructureFileReader reader, boolean alwaysConsiderOutgroupSampleGenotype) {
		this.svLocusGREventMap=new LinkedHashMap<>();
		
		for(SimpleSVLocus sv:this.coveredSVLocusList) {
			//find out all leaf nodes in subtree for each possible genotype; these leaf nodes on subtree may be non-leaf node on the original tree!!!!
			//notes that only Presence and absence genotypes can be found for leaf nodes in subtree
			Map<Genotype, List<TreeNode>> genotypeLeafNodesOnSubTreeMap = new HashMap<>();
			
			for(int id:this.tree.getNodeIDMap().keySet()) {
				TreeNode node = this.tree.getNodeIDMap().get(id);
				
				if(node.getSvLocusNodeInSubTreeMap().get(sv)) {//included in subtree
					if(node.isLeaf()) {//leaf node on original tree, always be included, trivial
						Genotype gt=node.getSvLocusGenotypeMap().get(sv);
						if(!genotypeLeafNodesOnSubTreeMap.containsKey(gt)) {
							genotypeLeafNodesOnSubTreeMap.put(gt, new ArrayList<>());
						}
						genotypeLeafNodesOnSubTreeMap.get(gt).add(node);
					}else {//non-leaf node on the original tree, check if it is leaf node on the subtree 
						//
						boolean childNodeIncluded=false;
						for(TreeNode child:node.getChildNodeList()) {
							if(child.getSvLocusNodeInSubTreeMap().get(sv)) {
								childNodeIncluded=true;
							}
						}
						if(!childNodeIncluded) {
							Genotype gt=node.getSvLocusGenotypeMap().get(sv);
							if(!genotypeLeafNodesOnSubTreeMap.containsKey(gt)) {
								genotypeLeafNodesOnSubTreeMap.put(gt, new ArrayList<>());
							}
							genotypeLeafNodesOnSubTreeMap.get(gt).add(node);
						}
					}
				}
			}
			
			
			if(alwaysConsiderOutgroupSampleGenotype) {
				/////////////////////for the two genotypes (0/0 and 1/1), 
				///////////////////////check if one genotype has 1 leaf node on subtree and the other genotype has 2 or more leaf node on tree
				int absenceGenotypeNum=genotypeLeafNodesOnSubTreeMap.containsKey(Genotype.ABSENCE)?genotypeLeafNodesOnSubTreeMap.get(Genotype.ABSENCE).size():0;
				int presenceGenotypeNum=genotypeLeafNodesOnSubTreeMap.containsKey(Genotype.PRESENCE)?genotypeLeafNodesOnSubTreeMap.get(Genotype.PRESENCE).size():0;
				
				if(absenceGenotypeNum==1 && presenceGenotypeNum>=2) {//GR event of the reverse of the SV occurred on the tree path starting from the node with absence genotype (0/0)
					//identify the tree branch path
					TreeNode start = genotypeLeafNodesOnSubTreeMap.get(Genotype.ABSENCE).get(0);
					List<TreeNode> treeNodesOnPath = new ArrayList<>();
					treeNodesOnPath.add(start);
					
					//trace back the path until the first ancestral node with one leaf node on sub-tree with genotype 'presence' (1/1) is encountered or the root node is reached
					TreeNode ancestralNode=start.getParent();
					if(ancestralNode!=null) {
						treeNodesOnPath.add(ancestralNode);
					}
					while(ancestralNode!=null && !ancestralNode.haveDescendantNodeInSubTreeWithGivenGenotypeForGivenSV(sv, Genotype.PRESENCE)) {//neither the ancestral node with the one leaf node on sub-tree with genotype 'presence' (1/1) is encountered nor the root node is reached
						ancestralNode=ancestralNode.getParent();
						if(ancestralNode!=null) {
							treeNodesOnPath.add(ancestralNode);
						}
					}
					
					//////////////////////////////////
					DescendantToAncestralTreePath treePath=new DescendantToAncestralTreePath(treeNodesOnPath);
					
					if(treePath.getBranchNumOnPath()>0) {//only add the GREvent if at least one branch is included
						GREvent grEvent=new GREvent(sv, this, true, treePath);
						
						if(DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromIngroupRootNodeToRootNode())
								||DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromOutgroupLeafNodeToRootNode())) {
							//skip if inferred GR event is only located on branches outside of the ingroup clade
						}else {
							this.svLocusGREventMap.put(sv, grEvent);
						}
					}
	
					
				}else if(absenceGenotypeNum>=2 && presenceGenotypeNum==1) {//GR event consistent with the SV occurred on the tree path starting from the node with presence genotype (0/0)
					//identify the tree branch path
					TreeNode start = genotypeLeafNodesOnSubTreeMap.get(Genotype.PRESENCE).get(0);
					List<TreeNode> treeNodesOnPath = new ArrayList<>();
					treeNodesOnPath.add(start);
					
					//trace back the path until the first ancestral node with one leaf node on sub-tree with genotype 'presence' (1/1) is encountered or the root node is reached
					TreeNode ancestralNode=start.getParent();
					if(ancestralNode!=null) {
						treeNodesOnPath.add(ancestralNode);
					}
					while(ancestralNode!=null && !ancestralNode.haveDescendantNodeInSubTreeWithGivenGenotypeForGivenSV(sv, Genotype.ABSENCE)) {//neither the ancestral node with the one leaf node on sub-tree with genotype 'presence' (1/1) is encountered nor the root node is reached
						ancestralNode=ancestralNode.getParent();
						if(ancestralNode!=null) {
							treeNodesOnPath.add(ancestralNode);
						}
					}
					//////////////////////////////////
					DescendantToAncestralTreePath treePath=new DescendantToAncestralTreePath(treeNodesOnPath);
					
					if(treePath.getBranchNumOnPath()>0) {//only add the GREvent if at least one branch is included
						GREvent grEvent=new GREvent(sv, this, false, treePath);
						
	
						if(DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromIngroupRootNodeToRootNode())
								||DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromOutgroupLeafNodeToRootNode())) {
							//skip if infered GR event is only located on branches outside of the ingroup clade
						}else {
							this.svLocusGREventMap.put(sv, grEvent);
						}
					}
				}else {//no GR event inferred
					//skip
				}
			}else {//only consider outgroup sample when necessary
				//find out the number of leaf node with with presence and absence genotype for all ingroup samples on the sub-tree
				//also find out the genotype of outgroup sample
				Genotype outgroupSampleGenotype=null;
				int ingroupSampleNumWithAbsenceGenotypeNum=0;
				int ingroupSampleNumWithPresenceGenotypeNum=0;
				
				if(genotypeLeafNodesOnSubTreeMap.containsKey(Genotype.ABSENCE)) {
					for(TreeNode node:genotypeLeafNodesOnSubTreeMap.get(Genotype.ABSENCE)) {
						//check if the tree node is within the clade containing all ingroup samples
						if(node.getParent()!=null) {//not root node
							//TODO
							if(reader.getAllOutgroupSampleIndices().contains(node.getSampleIndex())) {//outgroup sample node
								outgroupSampleGenotype=Genotype.ABSENCE;
							}else {//not root node, not outgroup sample node, must be within the clade of all ingroup samples
								ingroupSampleNumWithAbsenceGenotypeNum++;
							}
						}
					}
				}
				if(genotypeLeafNodesOnSubTreeMap.containsKey(Genotype.PRESENCE)) {
					for(TreeNode node:genotypeLeafNodesOnSubTreeMap.get(Genotype.PRESENCE)) {
						if(reader.getAllOutgroupSampleIndices().contains(node.getSampleIndex())) {//outgroup sample node
							outgroupSampleGenotype=Genotype.PRESENCE;
						}else {//not root node, not outgroup sample node, must be within the clade of all ingroup samples
							ingroupSampleNumWithPresenceGenotypeNum++;
						}
					}
				}
				/////////////////////////////////////
				if(ingroupSampleNumWithAbsenceGenotypeNum==1 && ingroupSampleNumWithPresenceGenotypeNum>=2 //
						|| (ingroupSampleNumWithAbsenceGenotypeNum==1 && ingroupSampleNumWithPresenceGenotypeNum==1)&&(outgroupSampleGenotype!=null && outgroupSampleGenotype.equals(Genotype.PRESENCE))
						////there is only one ingroup sample on the sub-tree for both presnce and absence genotype, check the outgroup sample genotype
						
						) {//GR event of the reverse of the SV occurred on the tree path starting from the node with absence genotype (0/0)
					//identify the tree branch path
					TreeNode start = genotypeLeafNodesOnSubTreeMap.get(Genotype.ABSENCE).get(0);
					List<TreeNode> treeNodesOnPath = new ArrayList<>();
					treeNodesOnPath.add(start);
					
					//trace back the path until the first ancestral node with one leaf node on sub-tree with genotype 'presence' (1/1) is encountered or the root node is reached
					TreeNode ancestralNode=start.getParent();
					if(ancestralNode!=null) {
						treeNodesOnPath.add(ancestralNode);
					}
					while(ancestralNode!=null && !ancestralNode.haveDescendantNodeInSubTreeWithGivenGenotypeForGivenSV(sv, Genotype.PRESENCE)) {//neither the ancestral node with the one leaf node on sub-tree with genotype 'presence' (1/1) is encountered nor the root node is reached
						ancestralNode=ancestralNode.getParent();
						if(ancestralNode!=null) {
							treeNodesOnPath.add(ancestralNode);
						}
					}
					//////////////////////////////////
					DescendantToAncestralTreePath treePath=new DescendantToAncestralTreePath(treeNodesOnPath);
					
					if(treePath.getBranchNumOnPath()>0) {//only add the GREvent if at least one branch is included
						GREvent grEvent=new GREvent(sv, this, true, treePath);
						
						if(DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromIngroupRootNodeToRootNode())
								||DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromOutgroupLeafNodeToRootNode())) {
							//skip if infered GR event is only located on branches outside of the ingroup clade
						}else {
							this.svLocusGREventMap.put(sv, grEvent);
						}
					}
					
				}else if(ingroupSampleNumWithAbsenceGenotypeNum>=2 && ingroupSampleNumWithPresenceGenotypeNum==1
						|| (ingroupSampleNumWithAbsenceGenotypeNum==1 && ingroupSampleNumWithPresenceGenotypeNum==1)&&(outgroupSampleGenotype!=null && outgroupSampleGenotype.equals(Genotype.ABSENCE))
						
						) {//GR event consistent with the SV occurred on the tree path starting from the node with presence genotype (0/0)
					//identify the tree branch path
					TreeNode start = genotypeLeafNodesOnSubTreeMap.get(Genotype.PRESENCE).get(0);
					List<TreeNode> treeNodesOnPath = new ArrayList<>();
					treeNodesOnPath.add(start);
					
					//trace back the path until the first ancestral node with one leaf node on sub-tree with genotype 'presence' (1/1) is encountered or the root node is reached
					TreeNode ancestralNode=start.getParent();
					if(ancestralNode!=null) {
						treeNodesOnPath.add(ancestralNode);
					}
					while(ancestralNode!=null && !ancestralNode.haveDescendantNodeInSubTreeWithGivenGenotypeForGivenSV(sv, Genotype.ABSENCE)) {//neither the ancestral node with the one leaf node on sub-tree with genotype 'presence' (1/1) is encountered nor the root node is reached
						ancestralNode=ancestralNode.getParent();
						if(ancestralNode!=null) {
							treeNodesOnPath.add(ancestralNode);
						}
					}

					//////////////////////////////////
					DescendantToAncestralTreePath treePath=new DescendantToAncestralTreePath(treeNodesOnPath);
					
					if(treePath.getBranchNumOnPath()>0) {//only add the GREvent if at least one branch is included
						GREvent grEvent=new GREvent(sv, this, false, treePath);
						
	
						if(DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromIngroupRootNodeToRootNode())
								||DescendantToAncestralTreePathUtils.sameTreePath(grEvent.getAssignedTreePath(), tree.getPathFromOutgroupLeafNodeToRootNode())) {
							//skip if infered GR event is only located on branches outside of the ingroup clade
						}else {
							this.svLocusGREventMap.put(sv, grEvent);
						}
					}
				}else {
					//skip
				}
			}
		}
	}
	
	//////////////////////////////////////////////
	
	
	/**
	 * @return the chrom
	 */
	public String getChrom() {
		return chrom;
	}


	/**
	 * @return the iD
	 */
	public int getID() {
		return ID;
	}

	/**
	 * @return the windowSize
	 */
	public int getWindowSize() {
		return windowSize;
	}

	/**
	 * @return the start
	 */
	public int getStart() {
		return start;
	}


	/**
	 * @return the end
	 */
	public int getEnd() {
		return end;
	}


	/**
	 * @return the tree
	 */
	public Tree getTree() {
		return tree;
	}

	/**
	 * @return the svLocusGREventMap
	 */
	public Map<SimpleSVLocus, GREvent> getSvLocusGREventMap() {
		return svLocusGREventMap;
	}

	/**
	 * @return the coveredSVLocusList
	 */
	public List<SimpleSVLocus> getCoveredSVLocusList() {
		return coveredSVLocusList;
	}
	
	
	
}
