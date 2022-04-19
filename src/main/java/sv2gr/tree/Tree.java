package sv2gr.tree;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import phylo.tree.reader.NewickFileFormatType;
import phylo.tree.reader.SimpleNewickParserUtils;


/**
 * tree class to facilitate GR inference from SV
 * 
 * @author tanxu
 *
 */
public class Tree {
	static int LEAF_INDEX;
	///////////////////////
	private TreeNode rootNode;
	
	/////////////////////////////
	private int counter;
	
	/**
	 * 
	 */
	private double longestRootLeafPathLength;
	
	/**
	 * largest number of edges between root and all leaf node
	 */
	private int maxLeafNodeEdgeNumToRoot;
	
	/**
	 * single 
	 */
	private double unitLength;
	
	private Integer totalLeafNum;
	
	
	Tree(){
		this.counter=0;
	}
	
	Tree(TreeNode rootNode){
		this.rootNode=rootNode;
		this.counter=0;
	}
	
	
	/**
	 * build and return the {@link DescendantToAncestralTreePath} from the root node of all ingroup samples to the root node of this tree
	 * 
	 * it is assumed there is only one single branch contained in the path
	 * 
	 * 
	 * @return
	 */
	public DescendantToAncestralTreePath getPathFromIngroupRootNodeToRootNode() {
		
		if(this.rootNode.getChildNodeList().size()!=2) {
			throw new IllegalArgumentException("root node does not contain two children nodes!");
		}
		
		List<TreeNode> treeNodesOnPath=new ArrayList<>();
		
		if(this.rootNode.getChildNodeList().get(0).isLeaf()) {//first child is outgroup leaf node
			treeNodesOnPath.add(this.rootNode.getChildNodeList().get(1));
		}else {
			treeNodesOnPath.add(this.rootNode.getChildNodeList().get(0));
		}
		
		treeNodesOnPath.add(this.rootNode);
		
		return new DescendantToAncestralTreePath(treeNodesOnPath);
	}
	
	
	/**
	 * build and return the {@link DescendantToAncestralTreePath} from the root node of all ingroup samples to the root node of this tree
	 * 
	 * it is assumed there is only one single branch contained in the path
	 * 
	 * @return
	 */
	public DescendantToAncestralTreePath getPathFromOutgroupLeafNodeToRootNode() {
		
		if(this.rootNode.getChildNodeList().size()!=2) {
			throw new IllegalArgumentException("root node does not contain two children nodes!");
		}
		
		List<TreeNode> treeNodesOnPath=new ArrayList<>();
		
		if(this.rootNode.getChildNodeList().get(0).isLeaf()) {//first child is outgroup leaf node
			treeNodesOnPath.add(this.rootNode.getChildNodeList().get(0));
		}else {
			treeNodesOnPath.add(this.rootNode.getChildNodeList().get(1));
		}
		
		treeNodesOnPath.add(this.rootNode);
		
		return new DescendantToAncestralTreePath(treeNodesOnPath);
	}
	
	
	//////////////////////////////////
	/**
	 * 
	 * @return
	 */
	int nextID() {
		int ret = counter;
		counter++;
		return ret;
	}
	
	void setRootNode(TreeNode rootNode) {
		this.rootNode = rootNode;
	}
	
	/**
	 * process the tree for visualization
	 * @param sortTreeCladesByIncreasingOrderOfLeafNum
	 */
	public void process(boolean sortTreeCladesByIncreasingOrderOfLeafNum) {
		//calculate tree and nodes attributes
		this.reorderChildrenNodes(sortTreeCladesByIncreasingOrderOfLeafNum);
		this.setAllLeafIndex();
		this.calculateNodeEdgeNum();
		this.calculateLargestDistBetweenRootAndLeaf();
	}
	
	
	/**
	 * reorder all children nodes of all nodes so that the one with less leaf nodes always come first
	 */
	public void reorderChildrenNodes(boolean increasingOrder) {
		this.rootNode.reorderChildrenNodesByCladeSize(increasingOrder);
	}
	
	
	/**
	 * 
	 */
	void setAllLeafIndex() {
		LEAF_INDEX=1;
		
		this.rootNode.setAllDescendantLeafIndex();
	}
	
	/**
	 * calculate the {@link #maxLeafNodeEdgeNumToRoot} and {@link #unitLength}
	 * 
	 * this will also first recursively invoke the {@link TreeNode#getMaxLeafNodeEdgeNumToThis()} method for all nodes on the tree from root to all leaves
	 */
	void calculateNodeEdgeNum() {
		
		this.maxLeafNodeEdgeNumToRoot=this.rootNode.getMaxLeafNodeEdgeNumToThis();
		
		this.unitLength = (double)1/this.maxLeafNodeEdgeNumToRoot;
	}
	
	/**
	 * calculate the {@link #longestRootLeafPathLength}
	 * 
	 * this will also invoke {@link TreeNode#getDistToRoot()} method for all nodes on the tree
	 */
	void calculateLargestDistBetweenRootAndLeaf() {
		this.longestRootLeafPathLength=0;
		
		for(int index:this.getNodeIDMap().keySet()) {
			if(this.getNodeIDMap().get(index).getDistToRoot()>this.longestRootLeafPathLength) {
				this.longestRootLeafPathLength=this.getNodeIDMap().get(index).getDistToRoot();
			}
		}
	}
	
	
	public int getTotalLeafNum() {
		if(this.totalLeafNum==null) {
			this.totalLeafNum=0;
			for(int index:this.getNodeIDMap().keySet()) {
				if(this.getNodeIDMap().get(index).isLeaf()) {
					this.totalLeafNum++;
				}
			}
		}
		return this.totalLeafNum;
	}

	
	
	
	////////////////////////////////
	/**
	 * @return the maxLeafNodeEdgeNumToRoot
	 */
	public int getMaxLeafNodeEdgeNumToRoot() {
		return maxLeafNodeEdgeNumToRoot;
	}
	
	/**
	 * @return the unitLength
	 */
	public double getUnitLength() {
		return unitLength;
	}


	/**
	 * @return the longestRootLeafPathLength
	 */
	public double getLongestRootLeafPathLength() {
		return longestRootLeafPathLength;
	}

	//////////////////////////////////
	
	/**
	 * build a tree with a single virtual root node from a newick tree string;
	 * 
	 * the given newickTreeStringFromFile must be in format '(...);' 
	 * 1. end with a ';'
	 * 2. enclosed by a pair of parenthesis that represent the virtual root node;
	 *  
	 * if not, a outer pair of parenthesis will be added before it is parsed to build the Tree object;
	 * 
	 * @param newickTreeStringFromFile the full newick tree string from a newick tree file, including the ending colon
	 * @param type
	 * @return
	 */
	public static Tree fromNewickString(String newickTreeStringFromFile, NewickFileFormatType type) {
		Tree ret = new Tree();
		TreeNode root = TreeNode.fromNewickTreeString(ret, SimpleNewickParserUtils.preprocessNewickStringFromFile(newickTreeStringFromFile), type, null);
		ret.setRootNode(root);
		return ret;
	}
	
	
	public static Tree fromNewickFile(Path singleTreeNewickFile, NewickFileFormatType type) {
		
		try {
			BufferedReader lineReader = new BufferedReader(new FileReader(singleTreeNewickFile.toFile()));
			String lineText = null;
			String newickTreeString="";
			while ((lineText = lineReader.readLine()) != null) {
				newickTreeString=newickTreeString.concat(lineText);
			}
			
			lineReader.close();
			
			return fromNewickString(newickTreeString, type);
		} catch (IOException ex) {
			ex.printStackTrace();
			return null;
		}
	}
	
	boolean isBifurcating() {
		return this.getRootNode().isBifurcating();
	}
	
	
	/**
	 * @return the rootNode
	 */
	public TreeNode getRootNode() {
		return rootNode;
	}
	
	/**
	 * return the nodes of this tree
	 * @return
	 */
	public Map<Integer, TreeNode> getNodeIDMap(){
		Map<Integer, TreeNode> ret = new LinkedHashMap<>();
		
		ret.put(this.rootNode.getId(), rootNode);
		
		ret.putAll(this.rootNode.getDescendantNodeIDMap());
		
		return ret;
	}
	
	/**
	 * 
	 * @return
	 */
	public Map<String, TreeNode> getLeafLabelTreeNodeMap(){
		Map<String, TreeNode> ret = new HashMap<>();
		
		Map<Integer, TreeNode> treeNodeIDMap = this.getNodeIDMap();
		for(int id:treeNodeIDMap.keySet()) {
			if(treeNodeIDMap.get(id).isLeaf()) {
				if(ret.containsKey(treeNodeIDMap.get(id).getLabel())) {
					throw new IllegalArgumentException("duplicate leaf label is found:"+treeNodeIDMap.get(id).getLabel());
				}
				ret.put(treeNodeIDMap.get(id).getLabel(), treeNodeIDMap.get(id));
			}
		}
		
		return ret;
	}
	
	/**
	 * return the set of nodes with average distance to all of its leaf nodes less than the given value
	 * 
	 * note that the return set of nodes should be disjoint with each other (every node is not the ancestor or descendant node of any other node)
	 * @param dist
	 * @return
	 */
	public Set<TreeNode> getTreeNodesWithAverageDistToAllLeafNodesLessThan(double dist){
		return this.rootNode.getDescendantNodesWithAverageDistToAllLeafNodesLessThan(dist);
	}
}
