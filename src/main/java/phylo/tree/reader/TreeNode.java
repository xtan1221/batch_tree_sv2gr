package phylo.tree.reader;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import basic.Pair;
import basic.Triple;


public class TreeNode {
	private final Tree tree;
	
	/**
	 * a unique integer among all nodes of the same tree;
	 */
	private final int id;
	private final TreeNode parent;
	
	////////////////
	
	///////////////
	private List<TreeNode> childNodeList;
	
	private TreeBranch branchToParent;
	
	/**
	 * index of this node
	 */
	private Integer leafIndex;
	
	////////////
	private Double distToParent;
	private Double bootstrap;
	private String label;
	
	
	/////////////////////////////////////
	/**
	 * number of edges from root node
	 */
	private Integer edgeNumToRoot;
	
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
	 * @return the branchToParent
	 */
	public TreeBranch getBranchToParent() {
		if(this.branchToParent==null) {
			if(this.parent!=null) {
				this.branchToParent=new TreeBranch(this.parent, this);
			}
		}
		return branchToParent;
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
        }else{//this node is an internal node with children nodes
            List<String> childrenNodeStringList = SimpleNewickParserUtils.splitNakedInternalNodeStringIntoChildrenNodeStrings(childrenNodesString);
            
            for(String childNodeString:childrenNodeStringList){
            	childNodeList.add(fromNewickTreeString(tree, childNodeString, type, ret));
            }
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
	
	@Override
	public String toString() {
		return "TreeNode [id=" + id + ", distToParent=" + distToParent + ", bootstrap=" + bootstrap + ", label=" + label
				+ "]";
	}

	/**
	 * @return the tree
	 */
	public Tree getTree() {
		return tree;
	}


	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + id;
		result = prime * result + ((parent == null) ? 0 : parent.hashCode());
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
		if (parent == null) {
			if (other.parent != null)
				return false;
		} else if (!parent.equals(other.parent))
			return false;
		return true;
	}
	
}
