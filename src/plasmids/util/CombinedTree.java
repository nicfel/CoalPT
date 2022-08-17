package plasmids.util;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.core.Input;
import beast.core.StateNode;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

public class CombinedTree extends StateNode implements TreeInterface {
	
    final public Input<List<Tree>> treeInput = new Input<>("tree", "trees that are supposed to be combined into one", new ArrayList<>());
    
    final public Input<List<TraitSet>> m_traitList = new Input<>("trait",
            "trait information for initializing traits (like node dates) in the tree",
            new ArrayList<>());
    final public Input<TaxonSet> m_taxonset = new Input<>("taxonset",
            "set of taxa that correspond to the leafs in the tree");
    final public Input<String> nodeTypeInput = new Input<>("nodetype",
            "type of the nodes in the beast.tree", Node.class.getName());
    
    /**
     * state of dirtiness of a node in the tree
     * DIRTY means a property on the node has changed, but not the topology. "property" includes the node height
     *       and that branch length to its parent.
     * FILTHY means the nodes' parent or child has changed.
     */
    public static final int IS_CLEAN = 0, IS_DIRTY = 1, IS_FILTHY = 2;



	@Override
	public void init(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public int getDimension() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getArrayValue(int dim) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public int getLeafNodeCount() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int getInternalNodeCount() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int getNodeCount() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public Node getRoot() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Node getNode(int i) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Node[] getNodesAsArray() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<Node> getExternalNodes() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<Node> getInternalNodes() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public TaxonSet getTaxonset() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void getMetaData(Node node, Double[] t, String pattern) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setMetaData(Node node, Double[] t, String pattern) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setEverythingDirty(boolean isDirty) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public StateNode copy() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void assignTo(StateNode other) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void assignFrom(StateNode other) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void assignFromFragile(StateNode other) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void fromXML(org.w3c.dom.Node node) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public int scale(double scale) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	protected void store() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void restore() {
		// TODO Auto-generated method stub
		
	}

}
