package conifer;

import java.util.List;

import org.jgrapht.UndirectedGraph;

import com.google.common.collect.Lists;

import fig.basic.UnorderedPair;



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
  
  public static <N> boolean isTip(UnorderedPair<N, N> edge, UndirectedGraph<N,UnorderedPair<N, N>> topology)
  {
    return topology.degreeOf(edge.getFirst()) == 1 || topology.degreeOf(edge.getSecond()) == 1;
  }
}
