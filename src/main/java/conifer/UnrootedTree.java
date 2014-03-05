package conifer;

import java.util.Map;

import org.jgraph.graph.DefaultEdge;
import org.jgrapht.UndirectedGraph;

import bayonet.graphs.GraphUtils;
import blang.annotations.Samplers;

import com.google.common.collect.Maps;

import conifer.moves.SingleBranchScaling;
import conifer.moves.SingleNNI;

import fig.basic.UnorderedPair;


@Samplers({SingleNNI.class, SingleBranchScaling.class})
public class UnrootedTree
{
  private final UndirectedGraph<TreeNode, UnorderedPair<TreeNode, TreeNode>> topology = GraphUtils.newUndirectedGraph();
  private final Map<UnorderedPair<TreeNode, TreeNode>,Double> branchLengths = Maps.newLinkedHashMap();
  
  public void addNode(TreeNode node)
  {
    this.topology.addVertex(node);
  }
  
  public void addEdge(TreeNode node1, TreeNode node2, double length)
  {
    this.topology.addEdge(node1, node2);
    branchLengths.put(new UnorderedPair<TreeNode, TreeNode>(node1, node2), length);
  }
  
  
  
  public UndirectedGraph<TreeNode, UnorderedPair<TreeNode, TreeNode>> getTopology()
  {
    return topology;
  }
  public Map<UnorderedPair<TreeNode, TreeNode>, Double> getBranchLengths()
  {
    return branchLengths;
  }
  
  public Double getBranchLength(TreeNode node1, TreeNode node2)
  {
    return branchLengths.get(new UnorderedPair<TreeNode, TreeNode>(node1, node2));
  }
  
  public void updateBranchLength(UnorderedPair<TreeNode, TreeNode> edge,
      double newValue)
  {
    if (!branchLengths.containsKey(edge))
      throw new RuntimeException();
    branchLengths.put(edge, newValue);
  }
  
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
