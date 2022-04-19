package sv2gr.gr;

import population.sv.utils.SimpleSVLocus;
import sv2gr.tree.DescendantToAncestralTreePath;
import sv2gr.tree.RegionalTree;
import sv2gr.tree.TreeNode;

/**
 * class for a genome rearrangement event
 * 
 * @author tanxu
 *
 */
public class GREvent {
	
	/////////////////////////////
	/**
	 * 
	 */
	private final SimpleSVLocus originalSVLocus;
	/**
	 * the RegionalTree on which this {@link GREvent} was inferred
	 */
	private final RegionalTree regionalTree;
	
	/**
	 * whether this GREvent is the reversal of the original SV
	 */
	private final boolean reverse;
	/**
	 * the path on the tree during which period this {@link GREvent} occurred
	 */
	private final DescendantToAncestralTreePath assignedTreePath;
	
	//////////////////////////////////
	/**
	 * number of branches on the {@link #assignedTreePath}
	 */
	private int branchNumOnDescendantToAncestralTreePath;
	/**
	 * the average distance from the ancestor node of the {@link #assignedTreePath} to all leaf nodes on the clade with the descendant node as root
	 * equal to the sum of {@link #descendantNodeAverageDistToLeafNodes} and length of {@link #assignedTreePath}
	 */
	private double ancestralNodeAverageDistToLeafNodes;
	
	/**
	 * the average distance from the descendant node of the {@link #assignedTreePath} to all leaf nodes on the clade with the descendant node as root
	 * equal to the {@link TreeNode#getAverageDistToAllDescendantLeafNodes()}
	 */
	private double descendantNodeAverageDistToLeafNodes;
	
	
	/**
	 * 
	 * @param originalSVLocus
	 * @param regionalTree
	 * @param reverse
	 * @param assignedTreePath
	 */
	public GREvent(SimpleSVLocus originalSVLocus, RegionalTree regionalTree, boolean reverse, DescendantToAncestralTreePath assignedTreePath) {
		super();
		this.originalSVLocus = originalSVLocus;
		this.regionalTree=regionalTree;
		this.reverse = reverse;
		this.assignedTreePath = assignedTreePath;
		
		/////////////////
		this.branchNumOnDescendantToAncestralTreePath=this.assignedTreePath.getBranchNumOnPath();
		
		this.descendantNodeAverageDistToLeafNodes=this.assignedTreePath.getDescendantNode().getAverageDistToAllDescendantLeafNodes();
		this.ancestralNodeAverageDistToLeafNodes = this.descendantNodeAverageDistToLeafNodes+this.assignedTreePath.getPathLen();
		
		//////////////////
		
	}
	

	/**
	 * @return the regionalTree
	 */
	public RegionalTree getRegionalTree() {
		return regionalTree;
	}


	/**
	 * @return the originalSVLocus
	 */
	public SimpleSVLocus getOriginalSVLocus() {
		return originalSVLocus;
	}

	
	/**
	 * @return the reverse
	 */
	public boolean isReversalOfSV() {
		return reverse;
	}

	
	/**
	 * @return the assignedTreePath
	 */
	public DescendantToAncestralTreePath getAssignedTreePath() {
		return assignedTreePath;
	}
	/////////////////////////


	/**
	 * @return the branchNumOnDescendantToAncestralTreePath
	 */
	public int getBranchNumOnDescendantToAncestralTreePath() {
		return branchNumOnDescendantToAncestralTreePath;
	}


	/**
	 * @return the ancestralNodeAverageDistToLeafNodes
	 */
	public double getAncestralNodeAverageDistToLeafNodes() {
		return ancestralNodeAverageDistToLeafNodes;
	}

	
	/**
	 * @return the descendantNodeAverageDistToLeafNodes
	 */
	public double getDescendantNodeAverageDistToLeafNodes() {
		return descendantNodeAverageDistToLeafNodes;
	}
	
	
}
