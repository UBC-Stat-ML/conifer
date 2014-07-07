package conifer;

import java.util.List;

import org.jgrapht.UndirectedGraph;

import bayonet.graphs.GraphUtils;
import briefj.BriefCollections;
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
      if (!GraphUtils.isTip(edge, topology))
        result.add(edge);
    
    return result;
  }
  
  public static TreeNode arbitraryNode(UnrootedTree tree)
  {
    return BriefCollections.pick(tree.getTopology().vertexSet());
  }
  
  public static List<TreeNode> syntheticTaxaList(int n)
  {
    List<TreeNode> result = Lists.newArrayList();
    for (int i = 0; i < n; i++)
      result.add(TreeNode.withLabel("synthetic-" + i));
    return result;
  }
  
}
