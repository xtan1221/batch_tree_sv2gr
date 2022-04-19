package phylo.tree.reader;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jgrapht.graph.SimpleGraph;

public class UnrootedTreeGraph {
	private SimpleGraph<GraphNode,GraphEdge> sug;
	
	private Map<Integer, GraphNode> treeNodeIdGraphNodeMap;
	
	/**
	 * build a undirected simple graph from a tree;
	 * 
	 * each node of the graph is a node on the tree;
	 * 
	 * the input tree must be bifurcating
	 * @param bifurnatingTree
	 */
	public UnrootedTreeGraph(Tree tree){
//		if(!bifurcatingTree.isBifurcating()) {
//			throw new IllegalArgumentException("given tree is not bifurcating!");
//		}
		
		this.sug = new SimpleGraph<>(GraphEdge.class);
		this.treeNodeIdGraphNodeMap = new HashMap<>();
		
		//for each node in the given tree (including root node), add a vertex to the graph 
		tree.getNodeIDMap().forEach((id, node)->{
//			if(node.getParent()!=null) {
			GraphNode graphNode = new GraphNode(id, node.getLabel());
			this.sug.addVertex(graphNode);
			treeNodeIdGraphNodeMap.put(id, graphNode);
//			}
		});
		
		////////////////////////////step 1
//		//add an edge between the two children nodes of the root node with the length equal to the sum of the distance to root node of the two children nodes;
//		TreeNode rootNodeChildNode1=tree.getRootNode().getChildNodeList().get(0);
//		TreeNode rootNodeChildNode2=tree.getRootNode().getChildNodeList().get(0);
//		
////		if(rootNodeChildNode1.getBootstrap()==null&&rootNodeChildNode2!=null || rootNodeChildNode1.getBootstrap()!=null && rootNodeChildNode2.getBootstrap()==null) {
////			throw new IllegalArgumentException("given tree's root node's two children node does not have consistent bootstrap!");
////		}else 
//		if(rootNodeChildNode1.getBootstrap()!=null && rootNodeChildNode2.getBootstrap()!=null && rootNodeChildNode1.getBootstrap()!=rootNodeChildNode2.getBootstrap()) {
//			throw new IllegalArgumentException("given tree's root node's two children node does not have equal bootstrap!");
//		}
//		
//		if(rootNodeChildNode1.getDistToParent()==null&&rootNodeChildNode2.getDistToParent()!=null || rootNodeChildNode1.getDistToParent()!=null && rootNodeChildNode2.getDistToParent()==null) {
//			throw new IllegalArgumentException("given tree's root node's two children node does not have consistent distance to parent!");
//		}
//		
//		
//		Double len = rootNodeChildNode1.getDistToParent()==null?null:rootNodeChildNode1.getDistToParent()+rootNodeChildNode2.getDistToParent();
//		Double bootstrap = null;
//		if(rootNodeChildNode1.getBootstrap()!=null) {
//			bootstrap=rootNodeChildNode1.getBootstrap();
//		}else if(rootNodeChildNode2.getBootstrap()!=null) {
//			bootstrap=rootNodeChildNode2.getBootstrap();
//		}
//		
//		GraphEdge graphEdge = new GraphEdge(
//				this.treeNodeIdGraphNodeMap.get(rootNodeChildNode1.getId()), 
//				this.treeNodeIdGraphNodeMap.get(rootNodeChildNode2.getId()), len, bootstrap);
//		
//		sug.addEdge(
//				this.treeNodeIdGraphNodeMap.get(rootNodeChildNode1.getId()), 
//				this.treeNodeIdGraphNodeMap.get(rootNodeChildNode2.getId()), 
//				graphEdge);
		
		///step 2. add other edges
		tree.getNodeIDMap().forEach((id, node)->{
			if(node.getParent()!=null) {
//				if(node.equals(rootNodeChildNode1)||node.equals(rootNodeChildNode2)) {
//					
//				}else {
					TreeNode parentNode = node.getParent();
					GraphEdge edge = new GraphEdge(
							this.treeNodeIdGraphNodeMap.get(node.getId()), 
							this.treeNodeIdGraphNodeMap.get(parentNode.getId()), 
							node.getDistToParent(), node.getBootstrap());
					sug.addEdge(
							this.treeNodeIdGraphNodeMap.get(node.getId()), 
							this.treeNodeIdGraphNodeMap.get(parentNode.getId()), 
							edge);
//				}
			}
		});
		
		
		
		///check if the tree's root node has exactly two children nodes, if yes, merge the two edges on the graph connecting the root node with its children nodes
		
		if(tree.getRootNode().getChildNodeList().size()==2) {
			GraphNode treeRootNodeOnGraph = this.treeNodeIdGraphNodeMap.get(tree.getRootNode().getId());
			assert this.sug.incomingEdgesOf(treeRootNodeOnGraph).size()==2;
			
			this.sug.removeVertex(treeRootNodeOnGraph);
			
			TreeNode childNode1=tree.getRootNode().getChildNodeList().get(0);
			TreeNode childNode2=tree.getRootNode().getChildNodeList().get(1);
			Double len = childNode1.getDistToParent()!=null?childNode1.getDistToParent()+childNode2.getDistToParent():null;
			Double bootstrap = childNode1.getBootstrap()!=null?childNode1.getBootstrap():childNode2.getBootstrap()!=null?childNode1.getBootstrap():null;
			
			GraphNode childNode1OnGraph=this.treeNodeIdGraphNodeMap.get(childNode1.getId());
			GraphNode childNode2OnGraph=this.treeNodeIdGraphNodeMap.get(childNode2.getId());
			
			GraphEdge newEdge = new GraphEdge(childNode1OnGraph, childNode2OnGraph, len, bootstrap);
			
			this.sug.addEdge(childNode1OnGraph, childNode2OnGraph, newEdge);
		}
	}
	
	
	/**
	 * build a Tree object by making the given leaf node as the outgroup
	 * 
	 * specifically, a new root node will be created in the middle of the outgroup node and its parent node in the original tree;
	 * 
	 * @param leaf
	 * @return
	 */
	Tree reroot(int outgroupLeafNodeID) {
		GraphNode outgroupLeaf = this.treeNodeIdGraphNodeMap.get(outgroupLeafNodeID);
		if(this.sug.incomingEdgesOf(outgroupLeaf).size()!=1) {
			throw new IllegalArgumentException("given node is not a leaf");
		}
		//////first modify the tree graph by adding the new root node and the edges
		GraphEdge adjacentEdge = sug.incomingEdgesOf(outgroupLeaf).iterator().next();
		
		assert adjacentEdge.node1.equals(outgroupLeaf)&&!adjacentEdge.node2.equals(outgroupLeaf) || !adjacentEdge.node1.equals(outgroupLeaf)&&adjacentEdge.node2.equals(outgroupLeaf);
		
		GraphNode adjacentNode = adjacentEdge.node1.equals(outgroupLeaf)?adjacentEdge.node2:adjacentEdge.node1;
		
		GraphNode newRootNode = new GraphNode(this.nextAvailableID(), null);
		
		this.sug.addVertex(newRootNode);
		this.sug.removeEdge(adjacentEdge);
		
		//edge between new root node and the outgroup leaf node; the bootstrap is always null
		GraphEdge edge1 = new GraphEdge(
				outgroupLeaf, newRootNode,
				adjacentEdge.len==null?null:adjacentEdge.len/2, 
				null//adjacentEdge.bootstrap==null?null:adjacentEdge.bootstrap
				);
		this.sug.addEdge(outgroupLeaf, newRootNode, edge1);
		
		//edge between new root node and the adjacentNode of the outgroup leaf node; the bootstrap is always null!
		GraphEdge edge2 = new GraphEdge(adjacentNode, newRootNode, adjacentEdge.len==null?null:adjacentEdge.len/2, null);
		this.sug.addEdge(adjacentNode, newRootNode, edge2);
		
		//////////
		//then 
		this.sug.vertexSet().forEach(v->{
			v.touched=false;
		});
		
		Tree ret = new Tree();
		TreeNode treeRootNode = new TreeNode(ret, newRootNode.id, null);
		ret.setRootNode(treeRootNode);
		
		treeRootNode.setBootstrap(null);
		treeRootNode.setDistToParent(null);
		treeRootNode.setLabel(newRootNode.label);
		newRootNode.touched=true;
		treeRootNode.setChildNodeList(this.buildChildrenTreeNodeList(ret, newRootNode, treeRootNode));
		
		
		return ret;
	}
	
	
	/**
	 * recursively build a list of TreeNode with the graph node adjacent to the given GraphNode and not touched yet;
	 * 
	 * @param graphNode
	 * @return
	 */
	private List<TreeNode> buildChildrenTreeNodeList(Tree tree, GraphNode graphNode, TreeNode treeNode){
		List<TreeNode> ret = new ArrayList<>();
		
		Set<GraphEdge> edges = this.sug.incomingEdgesOf(graphNode);
		
		//check every adjacent node of the given graphNode
		for(GraphEdge e: edges) {
			if(e.node1.touched&&e.node2.touched) {//node is already touched, skip
				//both node of the graph edge are touched, skip (the edge between the given node and its parent node on tree)
			}else {//node is not touched yet, must be a child tree node
				assert e.node1.touched&&!e.node2.touched || !e.node1.touched&&e.node2.touched;
				
				GraphNode childNode = e.node1.touched?e.node2:e.node1;
				
				TreeNode childTreeNode = new TreeNode(tree, childNode.id, treeNode);
				childTreeNode.setBootstrap(e.bootstrap);
				childTreeNode.setDistToParent(e.len);
				childTreeNode.setLabel(childNode.label);
				//must be set to true between 
				childNode.touched=true;
				
				//recursively build children nodes
				childTreeNode.setChildNodeList(this.buildChildrenTreeNodeList(tree, childNode, childTreeNode));
				
				ret.add(childTreeNode);
			}
		}
		
		return ret;
	}
	
	
	private int nextAvailableID() {
		int id=0;
		while(this.treeNodeIdGraphNodeMap.containsKey(id)) {
			id++;
		}
		return id;
	}
	
	
	static class GraphNode{
		private final int id;
		private final String label;
		private boolean touched = true;
		GraphNode(int id, String label){
			this.id = id;
			this.label=label;
		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + id;
			result = prime * result + ((label == null) ? 0 : label.hashCode());
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (!(obj instanceof GraphNode))
				return false;
			GraphNode other = (GraphNode) obj;
			if (id != other.id)
				return false;
			if (label == null) {
				if (other.label != null)
					return false;
			} else if (!label.equals(other.label))
				return false;
			return true;
		}
	}
	
	static class GraphEdge{
		private final GraphNode node1;
		private final GraphNode node2;
		private final Double len;
		private final Double bootstrap;
		GraphEdge(GraphNode node1, GraphNode node2, Double len, Double bootstrap){
			this.node1=node1;
			this.node2=node2;
			this.len=len;
			this.bootstrap=bootstrap;
		}
		
	}
}
