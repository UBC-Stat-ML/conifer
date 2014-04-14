package conifer;

import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.jgrapht.UndirectedGraph;

import bayonet.graphs.GraphUtils;
import bayonet.marginal.algo.EdgeSorter;
import blang.annotations.Processors;
import blang.annotations.Samplers;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;

import com.google.common.collect.Maps;

import conifer.moves.SingleBranchScaling;
import conifer.moves.SingleNNI;
import conifer.processors.TotalTreeLengthProcessor;
import conifer.processors.TreeDiameterProcessor;


/**
 * An unrooted phylogenetic tree.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
@Samplers({
  SingleNNI.class, 
  SingleBranchScaling.class
})
@Processors({
  TotalTreeLengthProcessor.class,
  TreeDiameterProcessor.class
})
public class UnrootedTree
{
  /**
   * Note if other fields are added, setTo() should be modified as well.
   */
  private UndirectedGraph<TreeNode, UnorderedPair<TreeNode, TreeNode>> topology = GraphUtils.newUndirectedGraph();
  private Map<UnorderedPair<TreeNode, TreeNode>,Double> branchLengths = Maps.newLinkedHashMap();
  
  /**
   * Add a node (species)
   * 
   * @param node
   */
  public void addNode(TreeNode node)
  {
    this.topology.addVertex(node);
  }
  
  /**
   * @param node1
   * @param node2
   * @param length
   */
  public void addEdge(TreeNode node1, TreeNode node2, double length)
  {
    this.topology.addEdge(node1, node2);
    branchLengths.put(new UnorderedPair<TreeNode, TreeNode>(node1, node2), length);
  }
  
  /**
   * 
   * @return
   */
  public UndirectedGraph<TreeNode, UnorderedPair<TreeNode, TreeNode>> getTopology()
  {
    return topology;
  }
  
  /**
   * 
   * @return
   */
  public Map<UnorderedPair<TreeNode, TreeNode>, Double> getBranchLengths()
  {
    return branchLengths;
  }
  
  /**
   * Return the branch corresponding to the unordered pair {node1, node2}
   * @param node1
   * @param node2
   * @return
   */
  public Double getBranchLength(TreeNode node1, TreeNode node2)
  {
    return branchLengths.get(new UnorderedPair<TreeNode, TreeNode>(node1, node2));
  }
  
  /**
   * Modifies in place the length of a branch.
   * @param edge
   * @param newValue
   */
  public void updateBranchLength(UnorderedPair<TreeNode, TreeNode> edge,
      double newValue)
  {
    if (!branchLengths.containsKey(edge))
      throw new RuntimeException();
    branchLengths.put(edge, newValue);
  }
  
  /**
   * Orient the edges relative to the provided root.
   * 
   * For each pair p, p.getLeft() contains the parent,
   * and p.getLeft() contains the children.
   * 
   * @param root
   * @return The oriented edges
   */
  public List<Pair<TreeNode,TreeNode>> getRootedEdges(TreeNode root)
  {
    return EdgeSorter.newEdgeSorter(getTopology(), root).backwardMessages();
  }
  
  /**
   * Performs a Nearest neighbor interchange move.
   * 
   * @param moved1
   * @param fixed1
   * @param moved2
   * @param fixed2
   */
  public void interchange(TreeNode moved1, TreeNode fixed1, TreeNode moved2, TreeNode fixed2)
  {
    modifyBranch(fixed1, moved1, fixed2);
    modifyBranch(fixed2, moved2, fixed1);
  }
  
  private void modifyBranch(
      TreeNode oldFixed, 
      TreeNode moved, 
      TreeNode newFixed)
  {
    topology.removeEdge(oldFixed, moved);
    topology.addEdge(moved, newFixed);
    double branchLength = getBranchLength(oldFixed, moved);
    branchLengths.remove(new UnorderedPair<TreeNode, TreeNode>(oldFixed, moved));
    branchLengths.put(new UnorderedPair<TreeNode,TreeNode>(moved, newFixed), branchLength);
  }

  /**
   * Set the value of the topology and branch length of this tree to those of the provided
   * tree (no copy performed, just two references set)
   * @param otherTree
   */
  public void setTo(UnrootedTree otherTree)
  {
    this.topology = otherTree.topology;
    this.branchLengths = otherTree.branchLengths;
  }

  /**
   * @see UnrootedTreeUtils.leaves()
   * @return
   */
  public List<TreeNode> leaves()
  {
    return TopologyUtils.leaves(getTopology());
  }
  
  /**
   * @see UnrootedTreeUtils.allTotalBranchLengthDistances()
   * @return
   */
  public Counter<UnorderedPair<TreeNode,TreeNode>> allTotalBranchLengthDistances()
  {
    return UnrootedTreeUtils.allTotalBranchLengthDistances(this);
  }
  
  public String toNewick()
  {
    return UnrootedTreeUtils.toNewick(this);
  }
  
  @Override
  public String toString()
  {
    return toNewick();
  }
}
