package sv2gr.tree;

import java.util.ArrayList;
import java.util.List;


/**
 * a path on tree from a node to one of its ancestral node (both ends inclusive)
 * 
 * @author tanxu
 */
public class DescendantToAncestralTreePath {
	private final List<TreeNode> nodesFromDescendantToAncestor;
	
	///////////////////
	/**
	 * length of this path; equal to the distance from the descendant node to the ancestor node
	 */
	private Double pathLength;
	/**
	 * the max value of bootstrap of all branches on this path
	 */
	private Double maxBoostrap;
	
	/**
	 * 
	 */
	public DescendantToAncestralTreePath(List<TreeNode> nodesFromDescendantToAncestor) {
		super();
		this.nodesFromDescendantToAncestor = nodesFromDescendantToAncestor;
	}

	/**
	 * @return the nodesFromDescendantToAncestor
	 */
	public List<TreeNode> getNodesFromDescendantToAncestor() {
		return nodesFromDescendantToAncestor;
	}

	/**
	 * 
	 * @return
	 */
	public TreeNode getDescendantNode() {
		return this.nodesFromDescendantToAncestor.get(0);
	}
	
	public TreeNode getAncestralNode() {
		return this.nodesFromDescendantToAncestor.get(this.nodesFromDescendantToAncestor.size()-1);
	}
	
	/**
	 * return whether this path does not contain any TreeNode
	 * @return
	 */
	public boolean isEmpty() {
		return this.nodesFromDescendantToAncestor.isEmpty();
	}
	
	
	/**
	 * return the number of branches on this path
	 * @return
	 */
	public int getBranchNumOnPath() {
		return this.nodesFromDescendantToAncestor.size()-1;
	}
	
	/**
	 * return the length of this path
	 * TODO testing
	 * @return
	 */
	public double getPathLen() {
		if(this.pathLength==null) {
			this.pathLength=0d;
			//add up the branch lengths
			for(int i=0;i<this.nodesFromDescendantToAncestor.size()-1;i++) {
				this.pathLength+=this.nodesFromDescendantToAncestor.get(i).getDistToParent();
			}
		}
		
		return this.pathLength;
	}
	
	/**
	 * 
	 * @return
	 */
	public Double getMaxBoostrap() {
		if(this.maxBoostrap==null) {
			for(int i=0; i<this.nodesFromDescendantToAncestor.size()-1;i++) {
				Double bs=this.nodesFromDescendantToAncestor.get(i).getBootstrap();
				if(this.maxBoostrap==null) {
					this.maxBoostrap=bs;
				}else {
					if(bs!=null && bs>this.maxBoostrap) {
						this.maxBoostrap=bs;
					}
				}
			}
		}
		return this.maxBoostrap;
	}
	
	/**
	 * return the string list of all leaf node sample names of the clade with the node next to the ancestral node as root
	 * 
	 * !!!!NOT with the ancestral node as root!!!
	 * @return
	 */
	public String getLeafSampleNameListStringOfAncestorNode() {
		StringBuilder sb=new StringBuilder();
		
		for(String sampleName:getLeafSampleNameListOfAncestorNode()) {
			sb.append(sampleName).append(":");
		}
		
		return sb.toString();
	}
	
	/**
	 * return the list of all leaf node sample names of the clade with the node next to the ancestral node as root
	 * 
	 * !!!!NOT with the ancestral node as root!!!
	 * 
	 * @return
	 */
	public List<String> getLeafSampleNameListOfAncestorNode() {
		List<String> ret = new ArrayList<>();
		
		TreeNode secondNodeNextToAncestralNode = this.nodesFromDescendantToAncestor.get(this.nodesFromDescendantToAncestor.size()-2);
		
		for(TreeNode leaf:secondNodeNextToAncestralNode.getAllDescendantLeafNodes()) {
			ret.add(leaf.getLabel());
		}
		
		return ret;
	}
	
	/**
	 * return the string list of all leaf node sample names of the clade with the descendant node as root
	 * 
	 * !!!!NOT with the ancestral node as root!!!
	 * @return
	 */
	public String getLeafSampleNameListStringOfDescendantNode() {
		StringBuilder sb=new StringBuilder();
				
		for(String sample:getLeafSampleNameListOfDescendantNode()) {
			sb.append(sample).append(":");
		}
		
		return sb.toString();
	}
	

	/**
	 * return the list of all leaf node sample names of the clade with the node next to the ancestral node as root
	 * 
	 * !!!!NOT with the ancestral node as root!!!
	 * 
	 * @return
	 */
	public List<String> getLeafSampleNameListOfDescendantNode() {
		List<String> ret = new ArrayList<>();
		
		TreeNode descendantNode = this.nodesFromDescendantToAncestor.get(0);
		
		for(TreeNode leaf:descendantNode.getAllDescendantLeafNodes()) {
			ret.add(leaf.getLabel());
		}
		
		return ret;
	}
}
