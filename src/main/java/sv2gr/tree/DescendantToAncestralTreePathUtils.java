package sv2gr.tree;

import java.util.ArrayList;
import java.util.List;

public class DescendantToAncestralTreePathUtils {
	
	
	
	/**
	 * check and return if the given two {@link DescendantToAncestralTreePath}s are on the same leaf-root path on the tree
	 * 
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static boolean onSameLeafRootPath(DescendantToAncestralTreePath p1, DescendantToAncestralTreePath p2) {
		if(p1.getDescendantNode().getAllAncestralNodes().containsAll(p2.getNodesFromDescendantToAncestor())
				||p2.getDescendantNode().getAllAncestralNodes().containsAll(p1.getNodesFromDescendantToAncestor())) {
			return true;
		}else {
			return false;
		}
	}
	
	
	/**
	 * check and return if p1 is fully located within p2
	 * 		note that this requires both of p1's boundary nodes should NOT be overlapping with p2's boundary tree nodes
	 * 
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static boolean p1IsFullyInsideOFP2(DescendantToAncestralTreePath p1, DescendantToAncestralTreePath p2) {
		if(!onSameLeafRootPath(p1,p2)) {
			return false;
		}
		
		if(p1.getDescendantNode().equals(p2.getDescendantNode())||p1.getAncestralNode().equals(p2.getAncestralNode())) {
			return false;
		}
		
		if(p2.getDescendantNode().getAllAncestralNodes().contains(p1.getDescendantNode()) 
				&& p1.getAncestralNode().getAllAncestralNodes().contains(p2.getAncestralNode())
				) {
			return true;
		}
		
		return false;
	}
	
	/**
	 * return the {@link DescendantToAncestralTreePath} shared by the given two {@link DescendantToAncestralTreePath}s
	 * 
	 * if the given two paths are not on the same leaf-root path, return null;
	 * if the given two paths are on same leaf-root path but does not share any branch, return empty {@link DescendantToAncestralTreePath}
	 * 
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static DescendantToAncestralTreePath getSharedTreePath(DescendantToAncestralTreePath p1, DescendantToAncestralTreePath p2) {
		if(p1.isEmpty()||p2.isEmpty()) {
			return new DescendantToAncestralTreePath(new ArrayList<>());
		}
		//
		if(!onSameLeafRootPath(p1, p2)) {
			return null;
		}
		
		List<TreeNode> treeNodeOnSharedPath = new ArrayList<>();
		if(p1.getDescendantNode().getAllAncestralNodes().containsAll(p2.getNodesFromDescendantToAncestor())) {
			TreeNode currentNode=p1.getDescendantNode();
			while(currentNode!=null) {
				if(p1.getNodesFromDescendantToAncestor().contains(currentNode) && p2.getNodesFromDescendantToAncestor().contains(currentNode)) {
					treeNodeOnSharedPath.add(currentNode);
				}
				currentNode = currentNode.getParent();
			}
		}else if(p2.getDescendantNode().getAllAncestralNodes().containsAll(p1.getNodesFromDescendantToAncestor())) {
			TreeNode currentNode=p2.getDescendantNode();
			while(currentNode!=null) {
				if(p1.getNodesFromDescendantToAncestor().contains(currentNode) && p2.getNodesFromDescendantToAncestor().contains(currentNode)) {
					treeNodeOnSharedPath.add(currentNode);
				}
				currentNode = currentNode.getParent();
			}
		}
		
		return new DescendantToAncestralTreePath(treeNodeOnSharedPath);
	}
	
	/**
	 * return the shortest {@link DescendantToAncestralTreePath} covering the given two {@link DescendantToAncestralTreePath}s
	 * 
	 * return null if there is no such {@link DescendantToAncestralTreePath}
	 * 
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static DescendantToAncestralTreePath getShortestCoveringTreePath(DescendantToAncestralTreePath p1, DescendantToAncestralTreePath p2) {
		if(!onSameLeafRootPath(p1, p2)) {
			return null;
		}
		
		List<TreeNode> treeNodeOnPath = new ArrayList<>();
		boolean p1AncestralNodeReached=false;
		boolean p2AncestralNodeReached=false;
		
		if(p1.getDescendantNode().getAllAncestralNodes().containsAll(p2.getNodesFromDescendantToAncestor())) {
			TreeNode currentNode=p1.getDescendantNode();
			treeNodeOnPath.add(currentNode);
			while(!(p1AncestralNodeReached&&p2AncestralNodeReached)) { //no need to check if currentNode is null or not, it will never be null before ancestral node of both paths are reached!
				currentNode = currentNode.getParent();
				treeNodeOnPath.add(currentNode);
				if(currentNode.equals(p1.getAncestralNode()))
					p1AncestralNodeReached=true;
				if(currentNode.equals(p2.getAncestralNode()))
					p2AncestralNodeReached=true;
			}
		}else if(p2.getDescendantNode().getAllAncestralNodes().containsAll(p1.getNodesFromDescendantToAncestor())) {
			TreeNode currentNode=p2.getDescendantNode();
			treeNodeOnPath.add(currentNode);
			while(!(p1AncestralNodeReached&&p2AncestralNodeReached)) {
				currentNode = currentNode.getParent();
				treeNodeOnPath.add(currentNode);
				if(currentNode.equals(p1.getAncestralNode()))
					p1AncestralNodeReached=true;
				if(currentNode.equals(p2.getAncestralNode()))
					p2AncestralNodeReached=true;
			}
		}
		
		return new DescendantToAncestralTreePath(treeNodeOnPath);
	}
	
	
	/**
	 * return whether the two given {@link DescendantToAncestralTreePath}s are the same
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static boolean sameTreePath(DescendantToAncestralTreePath p1, DescendantToAncestralTreePath p2) {
		return p1.getNodesFromDescendantToAncestor().equals(p2.getNodesFromDescendantToAncestor());
	}
	
	
}
