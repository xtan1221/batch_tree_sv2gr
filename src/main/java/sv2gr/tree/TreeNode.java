package sv2gr.tree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.base.Objects;

import basic.Pair;
import basic.Triple;
import phylo.tree.reader.NewickFileFormatType;
import phylo.tree.reader.SimpleNewickParserUtils;
import population.sv.utils.SimpleSVLocus;
import population.utils.Genotype;


public class TreeNode {
	private final Tree tree;
	
	/**
	 * a unique integer among all nodes of the same tree;
	 */
	private final int id;
	private final TreeNode parent;
	

	//////////////////////////////////
	/**
	 * index of this node
	 */
	private Integer leafIndex;
	
	////////////
	private Double distToParent;
	/**
	 * the percentage of bootstrap value in [0,1] if not null
	 * for leaf node, it is always assigned 1!!!
	 */
	private Double bootstrap; 
	private String label;
	
	
	///////////////
	private List<TreeNode> childNodeList;
	
	/**
	 * if this node is leaf, sample index contained by this leaf node
	 * otherwise, null
	 */
	private Integer sampleIndex;
	
	/**
	 * all covered {@link SimpleSVLocus} of the regional tree of this node
	 */
	private List<SimpleSVLocus> coveredSVLocusList;
	
	
	/**
	 * only contains the {@link SimpleSVLocus} with all ingroup samples either with 0/0 or 1/1 or ./. genotype
	 * cannot contain any {@link SimpleSVLocus} with at least one ingroup sample with 0/1 genotype
	 */
	private Map<SimpleSVLocus, Genotype> svLocusGenotypeMap;
	
	/**
	 * map from the SimpleSVLocus to whether this tree node is on the subtree of the SV for GR event inference
	 */
	private Map<SimpleSVLocus, Boolean> svLocusNodeInSubTreeMap;
	
	
	////////////////
	TreeNode(Tree tree, TreeNode parent){
		this.tree = tree;
		this.parent = parent;
		this.id = tree.nextID();
	}
	

	protected TreeNode(Tree tree, int index, TreeNode parent){
		this.tree = tree;
		this.parent = parent;
		this.id = index;
	}
	
	
	/**
	 * set the covered SVs of the regional tree of this node
	 * @param coveredSVLocusList the coveredSVLocusList to set
	 */
	public void setCoveredSVLocusList(List<SimpleSVLocus> coveredSVLocusList) {
		this.coveredSVLocusList = coveredSVLocusList;
	}
	
	/**
	 * set the sample index of this leaf node;
	 * 
	 * @param index
	 * @exception throw {@link UnsupportedOperationException} if this node is not leaf
	 */
	public void setSampleIndex(int index) {
		if(!this.isLeaf()) {
			throw new UnsupportedOperationException("tree node is not leaf!!!");
		}
		this.sampleIndex=index;
	}
	
	
	
//	/**
//	 * @return the ingroupSampleLeaf
//	 */
//	public Boolean getIngroupSampleLeaf() {
//		return ingroupSampleLeaf;
//	}
//
//	/**
//	 * @param ingroupSampleLeaf the ingroupSampleLeaf to set
//	 */
//	public void setIngroupSampleLeaf(Boolean ingroupSampleLeaf) {
//		this.ingroupSampleLeaf = ingroupSampleLeaf;
//	}
	
	
	/**
	 * @return the sampleIndex
	 */
	public Integer getSampleIndex() {
		return sampleIndex;
	}


	/**
	 * @return the svLocusGenotypeMap
	 */
	public Map<SimpleSVLocus, Genotype> getSvLocusGenotypeMap() {
		return svLocusGenotypeMap;
	}
	
	
	/**
	 * recursively infer the genotype of each SV on this TreeNode
	 * 
	 * if leaf node, directly assign the genotype of the corresponding sample from each SV
	 * 
	 * if internal nodes, 
	 * 		1. if all children nodes have the same non-null genotype (0/0,./.,1/1), set it to this node
	 * 			null genotype indicate that 
	 * 		2. else, set it to null
	 * 			either one or more child nodes have null genotype or they do not have the same non-null genotype
	 * 
	 * //TODO testing
	 */
	public void inferSVGenotypes() {
		//
		this.svLocusGenotypeMap = new HashMap<>();
		
		if(this.isLeaf()) {
			for(SimpleSVLocus sv:this.coveredSVLocusList) {
				this.svLocusGenotypeMap.put(sv, sv.getSampleIndexGenotypeMap().get(this.sampleIndex));
			}
		}else {
			//first invoke every child node's method
			for(TreeNode child:this.childNodeList) {
				child.inferSVGenotypes();
			}
			
			//assign genotype for each SV
			for(SimpleSVLocus sv:this.coveredSVLocusList) {
				//check if all child node have the same genotype (0/0 or 1/1, or ./.)
				Set<Genotype> childrenNodeNonNullGenotypes = new HashSet<>();
				boolean childNodeContainsNullGenotype = false;
				for(TreeNode child:this.childNodeList) {
					Genotype gt=child.getSvLocusGenotypeMap().get(sv);
					if(gt==null) {
						childNodeContainsNullGenotype = true;
					}else {
						childrenNodeNonNullGenotypes.add(gt);
					}
				}
				
				if(childNodeContainsNullGenotype || childrenNodeNonNullGenotypes.size()>=2) {//at least one child node has null genotype or all child nodes have multiple genotypes for the current sv 
					this.svLocusGenotypeMap.put(sv, null);
				}else {//a single non-null genotype
					this.svLocusGenotypeMap.put(sv, childrenNodeNonNullGenotypes.iterator().next());
				}
			}
		}
	}
	
	
	/**
	 * for each SV, check if this node should be included in the subtree for inference of GR events for that SV
	 * 
	 * 1. if a node is assigned non-null and non-missing genotype (0/0, 1/1)
	 * 		if its parent has the same genotype
	 * 			this node is NOT included 
	 * 		else 
	 * 			this node is included
	 * 2. if a node's genotype is null, it is included
	 * 
	 * 3. if a node's genotype is missing (./.), it is NOT included
	 * 
	 */
	public void inferInSubTree() {
		this.svLocusNodeInSubTreeMap=new HashMap<>();
		
		for(SimpleSVLocus sv:this.coveredSVLocusList) {
			Genotype gt=this.svLocusGenotypeMap.get(sv);
			if(gt==null) {
				this.svLocusNodeInSubTreeMap.put(sv, true);
			}else if(gt.equals(Genotype.MISSING)) {
				this.svLocusNodeInSubTreeMap.put(sv, false);
			}else {//non-null and non-missing, either 0/0 or 1/1
				Genotype parentGenotype = this.parent==null?null:this.parent.getSvLocusGenotypeMap().get(sv);
				
				if(Objects.equal(gt, parentGenotype)) {
					this.svLocusNodeInSubTreeMap.put(sv, false);
				}else {
					this.svLocusNodeInSubTreeMap.put(sv, true);
				}
			}
		}
		
		//invoke the child nodes' method
		if(!this.isLeaf()) {
			for(TreeNode child:this.childNodeList) {
				child.inferInSubTree();
			}
		}
	}
	

	
	/**
	 * @return the svLocusNodeInSubTreeMap
	 */
	public Map<SimpleSVLocus, Boolean> getSvLocusNodeInSubTreeMap() {
		return svLocusNodeInSubTreeMap;
	}

	///////////////////////
	

	//////////////////////////////////////////
//	/**
//	 * number of edges from root node
//	 */
//	private Integer edgeNumToRoot;
	
	private Double distToRoot;
	/**
	 * the average {@link #yPos} of all children nodes
	 * if this node is leaf, equal to the {@link #id}
	 */
	private Double yPos;
	
	/**
	 * the maximal edge number from descendant leaf node of this one to this node;
	 * for leaf node, this value is 0;
	 * 
	 */
	private Integer maxLeafNodeEdgeNumToThis;
	
	/**
	 * on a tree visualization with all real branch length ignored, all leaves have the same distance to root node
	 * 
	 * for a node, its x axis is Tree.
	 */
	private Double distToRootInBranchLengthIgnoredTree;
	
	////////////////
	/**
	 * all descendant leaf nodes of this one;
	 * if this is leaf node, itself is included in the {@link #allDescandentLeafNodes}
	 */
	private Set<TreeNode> allDescandentLeafNodes;
	
	/**
	 * the average distance from this node to all of its descendant leaf nodes
	 * 
	 * if this node is leaf, this value is 0;
	 */
	private Double averageDistToAllDescendantLeafNodes;
	
	
	
	/**
	 * reorder children nodes by number of their descendant leaf nodes
	 * this must be invoked before {@link #setAllDescendantLeafIndex} if reorder is needed
	 * @param increasingOrder whether sort children nodes by increasing order of number of leaf nodes
	 */
	public void reorderChildrenNodesByCladeSize(boolean increasingOrder) {
		if(this.isLeaf()) {//skip leaf nodes
			return;
		}
		
		Collections.sort(this.childNodeList, (a,b)->{
			if(increasingOrder) {
				return a.getAllDescendantLeafNodes().size()-b.getAllDescendantLeafNodes().size();
			}else {
				return b.getAllDescendantLeafNodes().size()-a.getAllDescendantLeafNodes().size();
			}
		});
		for(TreeNode child:this.childNodeList) {
			child.reorderChildrenNodesByCladeSize(increasingOrder);
		}
	}
	
	/**
	 * set {@link #leafIndex} of all descendant leaf nodes of this one
	 */
	public void setAllDescendantLeafIndex() {
		if(this.isLeaf()) {
			this.leafIndex = Tree.LEAF_INDEX;
			Tree.LEAF_INDEX++;
		}else {
			for(TreeNode child:this.childNodeList) {
				child.setAllDescendantLeafIndex();
			}
		}
	}
	
	/**
	 * recurively calcualte (if not already) and return {@link #maxLeafNodeEdgeNumToThis}
	 * 
	 * @return
	 */
	public int getMaxLeafNodeEdgeNumToThis() {
		if(maxLeafNodeEdgeNumToThis==null) {
			if(this.isLeaf()) {
				this.maxLeafNodeEdgeNumToThis=0;
			}else {
				int maxLeafNodeEdgeNumToChildNode=0;
				for(TreeNode child:this.childNodeList) {
					if(child.getMaxLeafNodeEdgeNumToThis()>maxLeafNodeEdgeNumToChildNode) {
						maxLeafNodeEdgeNumToChildNode=child.getMaxLeafNodeEdgeNumToThis();
					}
				}
				
				this.maxLeafNodeEdgeNumToThis=maxLeafNodeEdgeNumToChildNode+1;
			}
		}
		
		return this.maxLeafNodeEdgeNumToThis;
	}
	
	/**
	 * calculate and return the {@link #distToRootInBranchLengthIgnoredTree}
	 * can only be invoked after {@link #maxLeafNodeEdgeNumToThis} is calculated by {@link #getMaxLeafNodeEdgeNumToThis()}
	 * @return
	 */
	public double getDistToRootInBranchLengthIgnoredTree() {
		if(this.distToRootInBranchLengthIgnoredTree==null) {
			this.distToRootInBranchLengthIgnoredTree = this.tree.getUnitLength()*(this.tree.getMaxLeafNodeEdgeNumToRoot()-this.maxLeafNodeEdgeNumToThis);
		}
		return this.distToRootInBranchLengthIgnoredTree;
	}
	
	
	/**
	 * calculate and return the {@link #distToRoot}
	 * @return the distToRoot
	 */
	public Double getDistToRoot() {
		if(this.distToRoot==null) {
			if(this.getParent()==null) {
				this.distToRoot=this.distToParent;
			}else {
				this.distToRoot = this.distToParent+this.getParent().getDistToRoot();
			}
		}
		return this.distToRoot;
	}

	public double getYPos() {
		if(this.yPos==null) {
			if(this.isLeaf()) {
				return this.leafIndex;
			}else {
				double summedChildNodeYPos=0;
				for(TreeNode child:this.childNodeList) {
					summedChildNodeYPos+=child.getYPos();
				}
				this.yPos=summedChildNodeYPos/this.childNodeList.size();
			}
		}
		
		return yPos;
	}
	
	
	/**
	 * 
	 * @return
	 */
	public int getTotalDescandentLeafNodeNum() {
		return this.getAllDescendantLeafNodes().size();
	}
	
	
	public Set<TreeNode> getAllDescendantLeafNodes(){
		if(this.allDescandentLeafNodes==null) {
			this.allDescandentLeafNodes = new HashSet<>();
			if(this.isLeaf()) {
				this.allDescandentLeafNodes.add(this);
			}else {
				for(TreeNode child:this.childNodeList) {
					this.allDescandentLeafNodes.addAll(child.getAllDescendantLeafNodes());
				}
			}
		}
		
		return this.allDescandentLeafNodes;
	}
	
	/**
	 * all ancestral nodes of this one on the path from this node to the root node
	 */
	private List<TreeNode> ancestralNodes;
	
	/**
	 * 
	 * @return
	 */
	public List<TreeNode> getAllAncestralNodes(){
		if(this.ancestralNodes==null) {
			this.ancestralNodes=new ArrayList<>();
			
			if(this.parent==null) {
				
			}else {
				this.ancestralNodes.add(parent);
				this.ancestralNodes.addAll(parent.getAllAncestralNodes());
			}
		}
		
		return this.ancestralNodes;
	}
	
	/**
	 * 
	 * @return
	 */
	public double getAverageDistToAllDescendantLeafNodes() {
		if(this.averageDistToAllDescendantLeafNodes==null) {
			if(this.isLeaf()) {
				this.averageDistToAllDescendantLeafNodes=0d;
			}else {
				double summedChild = 0;
				
				for(TreeNode child:this.childNodeList) {
					summedChild=summedChild+(child.getAverageDistToAllDescendantLeafNodes()+child.getDistToParent())*child.getTotalDescandentLeafNodeNum();
				}
				
				this.averageDistToAllDescendantLeafNodes = summedChild/this.getTotalDescandentLeafNodeNum();
			}
		}
		
		return this.averageDistToAllDescendantLeafNodes;
	}
	
	
	
	/**
	 * find and return 
	 * @param dist
	 * @return
	 */
	public Set<TreeNode> getDescendantNodesWithAverageDistToAllLeafNodesLessThan(double dist){
		Set<TreeNode> ret = new HashSet<>();
		
		if(this.isLeaf()) {
			ret.add(this);
		}else {
			for(TreeNode child:this.childNodeList) {
				if(child.getAverageDistToAllDescendantLeafNodes()<dist) {
					ret.add(child);
				}else {
					ret.addAll(child.getDescendantNodesWithAverageDistToAllLeafNodesLessThan(dist));
				}
			}
		}
		
		return ret;
	}
	
	/**
	 * 
	 * @return
	 */
	public Set<TreeNode> getAllDescendantNodes(){
		Set<TreeNode> ret = new HashSet<>();
		
		if(this.isLeaf()) {
			ret.add(this);
		}else {
			for(TreeNode child:this.childNodeList) {
				ret.add(child);
				ret.addAll(child.getAllDescendantNodes());
			}
		}
		
		return ret;
	}
	
	/**
	 * check if this nodes have at least one descendant node that is in subtree and with the given genotype for the given SV
	 * 
	 * @param sv
	 * @param gt
	 * @return
	 */
	public boolean haveDescendantNodeInSubTreeWithGivenGenotypeForGivenSV(SimpleSVLocus sv, Genotype gt) {
		if(this.svLocusNodeInSubTreeMap.get(sv)) {//this node is in subtree for the given sv
			if(Objects.equal(this.svLocusGenotypeMap.get(sv), gt)) {//this node has the same genotype as the given one
				return true;
			}else {//check child node recursively
				for(TreeNode child:this.childNodeList) {
					if(child.haveDescendantNodeInSubTreeWithGivenGenotypeForGivenSV(sv, gt)) {
						return true;
					}
				}
				
				return false;
			}
		}else {//this node is not in subtree for the given sv
			return false;
		}
	}
	
	/////////////////////////////////////////////////////////
	
	/**
	 * 
	 * @param newickTreeString
	 * @param type
	 * @param parent
	 * @return
	 */
	static TreeNode fromNewickTreeString(Tree tree, String newickTreeString, NewickFileFormatType type, TreeNode parent) {
		TreeNode ret = new TreeNode(tree, parent);
		
		Triple<String,String,String> splits = SimpleNewickParserUtils.extractChildrenNodeStringNodeLabelStringAndBranchLabelString(newickTreeString, type);
        String childrenNodesString = splits.getLeft();
        
        String labelString = splits.getMiddle();
        ret.label = labelString;
        
        //
        String edgeLabelString = splits.getRight(); //for edge length and bootstrap
        Pair<Double,Double> parsedResult = SimpleNewickParserUtils.parseEdgeLabelStringForLengthAndBootstrap(edgeLabelString, type);
        Double distToParent = parsedResult.getFirst();
        Double boostrapToParent = parsedResult.getSecond();
        ret.distToParent=distToParent==null?0:distToParent;
        ret.bootstrap=boostrapToParent;
        
        //
        List<TreeNode> childNodeList = new ArrayList<>();
        if(childrenNodesString == null){//this node is a leaf
            //no children
        	ret.bootstrap=1d; //leaf node's branch to parent node is always 100% bootstrapped!!!!
        }else{//this node is an internal node with children nodes
            List<String> childrenNodeStringList = SimpleNewickParserUtils.splitNakedInternalNodeStringIntoChildrenNodeStrings(childrenNodesString);
            
            for(String childNodeString:childrenNodeStringList){
            	childNodeList.add(fromNewickTreeString(tree, childNodeString, type, ret));
            }
        }
        
        ///the branch from child node of root node to root node should always have bootstrap = 1!!!!!!!!!!
        if(ret.getParent()!=null && ret.getParent().getParent()==null) {//current node is a child node of the root node
        	ret.bootstrap=1d;
        }
        
		ret.childNodeList=childNodeList;
		
		
		
		return ret;
	}
	
	
	/**
	 * build and return a newick tree string of this Node
	 * =================
	 * SIMPLE_NEWICK_1("SIMPLE_NEWICK_1", new VfNotes("bootstrap value (if exist) is in squared bracket after branch length like in \n\t ((raccoon:19.19959,bear:6.80041):0.84600[50],((sea_lion:11.99700, seal:12.00300):7.52973[100],((monkey:100.85930,cat:47.14069):20.59201[80], weasel:18.87953):2.09460[75]):3.87382[50],dog:25.46154);")),  
	 * 	//this type of newick seems not allow internal node labels since the edge length takes the position of the internal node label of SIMPLE_NEWICK_1?
	 * 	//thus parsing of this type regarding the 
	 * 	SIMPLE_NEWICK_2("SIMPLE_NEWICK_2", new VfNotes("bootstrap value (if exist) is before branch length with a colon like in \n\t ((raccoon:19.19959,bear:6.80041)50:0.84600,((sea_lion:11.99700, seal:12.00300)100:7.52973,((monkey:100.85930,cat:47.14069)80:20.59201, weasel:18.87953)75:2.09460)50:3.87382,dog:25.46154);")),  
	 * @return
	 */
	String toNewickString(NewickFileFormatType type) {
		String ret = "";
		
		if(this.childNodeList.isEmpty()) {//leaf raccoon:19.19959 or raccoon, no parenthesis, no bootstrap
			ret = this.label;
			if(this.distToParent!=null)
				ret = ret.concat(":").concat(Double.toString(distToParent));
		}else {//internal node
			ret="(";
			boolean firstAdded = false;
			for(TreeNode child:this.childNodeList) {
				if(firstAdded) {
					ret = ret.concat(",");
				}else {
					firstAdded = true;
				}
				
				ret = ret.concat(child.toNewickString(type));
			}
			
			ret = ret.concat(")");
			
			if(this.bootstrap!=null) {
				if(type.equals(NewickFileFormatType.SIMPLE_NEWICK_1)) {//(raccoon:19.19959,bear:6.80041):0.84600[50] or (raccoon:19.19959,bear:6.80041)[50] 
					if(this.distToParent!=null) {//(raccoon:19.19959,bear:6.80041):0.84600[50]
						ret=ret.concat(":").concat(Double.toString(this.distToParent)).concat("[").concat(Double.toString(this.bootstrap)).concat("]");
					}else {//(raccoon:19.19959,bear:6.80041)[50] 
						ret=ret.concat("[").concat(Double.toString(this.bootstrap)).concat("]");
					}
				}else {//(raccoon:19.19959,bear:6.80041)50:0.84600 or (raccoon:19.19959,bear:6.80041)50
					if(this.distToParent!=null) {//(raccoon:19.19959,bear:6.80041)50:0.84600
						ret=ret.concat(Double.toString(this.bootstrap)).concat(":").concat(Double.toString(this.distToParent));
					}else {//(raccoon:19.19959,bear:6.80041)50
						ret=ret.concat(Double.toString(this.bootstrap));
					}
				}
			}else {//no bootstrap must be (raccoon:19.19959,bear:6.80041):0.84600 or (raccoon:19.19959,bear:6.80041)
				if(this.distToParent!=null) {//(raccoon:19.19959,bear:6.80041):0.84600
					ret=ret.concat(":").concat(Double.toString(this.distToParent));
				}else {//(raccoon:19.19959,bear:6.80041)
					//do nothing
				}
			}
		}
		
		return ret;
	}
	
	/**
	 * @return the id
	 */
	public int getId() {
		return id;
	}

	/**
	 * @return the parent
	 */
	public TreeNode getParent() {
		return parent;
	}

	/**
	 * @return the childNodeList
	 */
	public List<TreeNode> getChildNodeList() {
		return childNodeList;
	}

	/**
	 * @return the distToParent
	 */
	public Double getDistToParent() {
		return distToParent;
	}

	/**
	 * @return the bootstrap
	 */
	public Double getBootstrap() {
		return bootstrap;
	}

	/**
	 * @return the label
	 */
	public String getLabel() {
		return label;
	}
	
	

	/**
	 * @param childNodeList the childNodeList to set
	 */
	public void setChildNodeList(List<TreeNode> childNodeList) {
		this.childNodeList = childNodeList;
	}



	/**
	 * @param distToParent the distToParent to set
	 */
	public void setDistToParent(Double distToParent) {
		this.distToParent = distToParent;
	}



	/**
	 * @param bootstrap the bootstrap to set
	 */
	public void setBootstrap(Double bootstrap) {
		this.bootstrap = bootstrap;
	}



	/**
	 * @param label the label to set
	 */
	public void setLabel(String label) {
		this.label = label;
	}

	
	/**
	 * return the descendant nodes of this Node (exclusive)
	 * @return
	 */
	Map<Integer, TreeNode> getDescendantNodeIDMap(){
		Map<Integer, TreeNode> ret = new LinkedHashMap<>();
		
		this.childNodeList.forEach(c->{
			ret.put(c.getId(), c);
			ret.putAll(c.getDescendantNodeIDMap());
		});
		
		return ret;
	}
	
	public boolean isLeaf() {
		return this.childNodeList.size()==0;
	}
	
	boolean isBifurcating() {
		if(this.isLeaf()) {
			return true;
		}else if(this.getChildNodeList().size()==2) { //
			boolean ret = true;
			for(TreeNode child:this.getChildNodeList()) {
				ret = ret&&child.isBifurcating();
			}
			return ret;
		}else {//non-bifurcating
			return false;
		}
	}
	
	

	/**
	 * @return the tree
	 */
	public Tree getTree() {
		return tree;
	}
	
	
	
	/////////////////////////////
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + id;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!(obj instanceof TreeNode))
			return false;
		TreeNode other = (TreeNode) obj;
		if (id != other.id)
			return false;
		return true;
	}

	@Override
	public String toString() {
		return "TreeNode [id=" + id + ", distToParent=" + distToParent + ", bootstrap=" + bootstrap + ", label=" + label
				+ "]";
	}

}
