package conifer;

import java.util.List;

import org.jgrapht.UndirectedGraph;

import briefj.collections.UnorderedPair;

import com.google.common.collect.Lists;



/**
 * Utilities related to tree topologies.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class TopologyUtils
{
  /**
   * @param <N>
   * @param topology
   * @return Edges that are not tips
   */
  public static <N> List<UnorderedPair<N, N>> nonTerminalEdges(UndirectedGraph<N,UnorderedPair<N, N>> topology)
  {
    List<UnorderedPair<N, N>> result = Lists.newArrayList();
    
    for (UnorderedPair<N, N> edge : topology.edgeSet())
      if (!isTip(edge, topology))
        result.add(edge);
    
    return result;
  }
  
  /**
   * @param <N>
   * @param edge
   * @param topology
   * @return The edges connected to leaves.
   */
  public static <N> boolean isTip(UnorderedPair<N, N> edge, UndirectedGraph<N,UnorderedPair<N, N>> topology)
  {
    return topology.degreeOf(edge.getFirst()) == 1 || topology.degreeOf(edge.getSecond()) == 1;
  }
  
  /**
   * 
   * @param <N>
   * @param topology
   * @return The nodes that have strictly more than one neighbors.
   */
  public static <N> List<N> internalNodes(UndirectedGraph<N,UnorderedPair<N, N>> topology)
  {
    List<N> result = Lists.newArrayList();
    for (N vertex : topology.vertexSet())
      if (topology.degreeOf(vertex) > 1)
        result.add(vertex);
    return result;
  }
  
  /**
   * 
   * @param <N>
   * @param topology
   * @return The nodes that have one or zero neighbors.
   */
  public static <N> List<N> leaves(UndirectedGraph<N,UnorderedPair<N, N>> topology)
  {
    List<N> result = Lists.newArrayList();
    for (N vertex : topology.vertexSet())
      if (topology.degreeOf(vertex) <= 1)
        result.add(vertex);
    return result;
  }
  
  public static List<TreeNode> syntheticTaxaList(int n)
  {
    List<TreeNode> result = Lists.newArrayList();
    for (int i = 0; i < n; i++)
      result.add(TreeNode.withLabel("synthetic-" + i));
    return result;
  }
  
}
