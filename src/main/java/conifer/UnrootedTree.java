package conifer;

import java.util.Map;

import org.jgrapht.UndirectedGraph;

import bayonet.graphs.GraphUtils;
import blang.annotations.Samplers;
import briefj.collections.UnorderedPair;

import com.google.common.collect.Maps;

import conifer.moves.SingleBranchScaling;
import conifer.moves.SingleNNI;


/**
 * An unrooted phylogenetic tree.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
@Samplers({SingleNNI.class, SingleBranchScaling.class})
public class UnrootedTree
{
  private final UndirectedGraph<TreeNode, UnorderedPair<TreeNode, TreeNode>> topology = GraphUtils.newUndirectedGraph();
  private final Map<UnorderedPair<TreeNode, TreeNode>,Double> branchLengths = Maps.newLinkedHashMap();
  
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
   * Performs a Nearest neighbor interchange move.
   * 
   * TODO: add picture
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

  
  
}
