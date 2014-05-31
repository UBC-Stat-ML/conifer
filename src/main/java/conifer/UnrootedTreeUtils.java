package conifer;


import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.lang3.tuple.Pair;
import org.jgrapht.Graphs;

import com.google.common.collect.Lists;

import conifer.factors.NonClockTreePrior;
import conifer.io.newick.NewickParser;
import conifer.io.newick.ParseException;

import bayonet.distributions.Exponential;
import bayonet.graphs.GraphUtils;
import blang.variables.RealVariable;
import briefj.BriefCollections;
import briefj.BriefIO;
import briefj.BriefLists;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.Tree;
import briefj.collections.UnorderedPair;



public class UnrootedTreeUtils
{
  /**
   * See comments in UnrootedTree
   * @param f
   * @return
   */
  public static UnrootedTree fromNewick(File f)
  {
    return fromNewickString(BriefIO.fileToString(f));
  }
  
  /**
   * See comments in UnrootedTree
   * @param f
   * @return
   */
  public static UnrootedTree fromNewickString(String newickString)
  {
    NewickParser parser = new NewickParser(newickString); 
    try
    {
      Tree<TreeNode> topo = parser.parse();
      Map<TreeNode,Double> bls = parser.getBranchLengths();
      UnrootedTree result = new UnrootedTree();
      process(result, topo, bls, null);
      return result;
    } catch (ParseException e)
    {
      throw new RuntimeException(e);
    }
  }
  
  private static void process(
      UnrootedTree result, Tree<TreeNode> topo,
      Map<TreeNode, Double> bls, TreeNode ancestor)
  {
    TreeNode current = topo.getLabel();
    result.addNode(current);
    if (ancestor != null)
      result.addEdge(ancestor, current, bls.get(current));
    for (Tree<TreeNode> child : topo.getChildren())
      process(result, child, bls, current);
  }

  public static String toNewick(UnrootedTree tree)
  {
    StringBuilder result = new StringBuilder();
    TreeNode pseudoRoot = null;
    // root the tree at an internal node if possible
    List<TreeNode> internalNodes = TopologyUtils.internalNodes(tree.getTopology());
    if (!internalNodes.isEmpty())
      pseudoRoot = BriefCollections.pick(internalNodes);
    else if (!tree.getTopology().vertexSet().isEmpty()) // otherwise, do it arbitrarily
      pseudoRoot = BriefCollections.pick(tree.getTopology().vertexSet());
    if (pseudoRoot != null)
      toNewick(tree, null, pseudoRoot, result);
    result.append(";");
    return result.toString();
  }
  
  private static void toNewick(UnrootedTree tree, TreeNode parent, TreeNode current, StringBuilder builder)
  {
    List<TreeNode> children = Graphs.neighborListOf(tree.getTopology(), current);
    if (parent != null)
      if (!children.remove(parent))
        throw new RuntimeException();
    
    if (!children.isEmpty())
    {
      builder.append("(");
      for (int cIndex = 0; cIndex < children.size(); cIndex++)
      {
        toNewick(tree, current, children.get(cIndex), builder);
        if (cIndex != children.size() - 1)
          builder.append(",");
      }
      builder.append(")");
    }
    if (current.isLabelled())
    {
      String label = current.toString();
      if (label.contains("(") || 
          label.contains(")") || 
          label.contains(",") || 
          label.contains(":") ||
          label.contains(";"))
        throw new RuntimeException();
      builder.append(label);
    }
    if (parent != null)
      builder.append(":" + tree.getBranchLength(current, parent));
  }
  
  /**
   * Compute the total branch length distance between each pair of leaves in a tree.
   * 
   * By total branch length distance, we mean the sum of branch lengths encountered
   * in the unique shortest path between two leaves.
   * 
   * @param tree
   * @return
   */
  public static Counter<UnorderedPair<TreeNode,TreeNode>> allTotalBranchLengthDistances(UnrootedTree tree)
  {
    return new EfficientUnrootedTree(tree).allTotalBranchLengthDistances();
  }
  
  /**
   * An efficient array-based implementation of unrooted trees used
   * internally to compute pairwise distances.
   * 
   * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
   *
   */
  private static class EfficientUnrootedTree
  {
    private final int [][] nbhrs;
    private final double [][] bls;
    private final int nLeaves;
    private final Indexer<TreeNode> indexer = new Indexer<TreeNode>();
    
    public EfficientUnrootedTree(UnrootedTree ut)
    {
      final List<TreeNode> leaves = ut.leaves();
      this.nLeaves = leaves.size();
      for (final TreeNode leaf : leaves)
        indexer.addToIndex(leaf);
      for (final TreeNode node : ut.getTopology().vertexSet())
        if (!indexer.containsObject(node))
          indexer.addToIndex(node);
      final int size = indexer.size();
      this.nbhrs = new int[size][];
      this.bls = new double[size][];
      for (final TreeNode t1 : ut.getTopology().vertexSet())
      {
        final int i1 = indexer.o2i(t1);
        int currentDegree = ut.getTopology().degreeOf(t1);
        nbhrs[i1] = new int[currentDegree];
        bls[i1] = new double[currentDegree];
        int idx = 0;
        for (TreeNode t2 : Graphs.neighborListOf(ut.getTopology(), t1))
        {
          final int i2 = indexer.o2i(t2);
          nbhrs[i1][idx] = i2;
          bls[i1][idx] = ut.getBranchLength(t1, t2);
          idx++;
        }
      }
    }
    
    public Counter<UnorderedPair<TreeNode,TreeNode>> allTotalBranchLengthDistances()
    {
      double [][] result = new double[nLeaves][nLeaves];
      
      for (int start = 0; start < nLeaves - 1; start++)
        _dfsTotalBL(0.0, start, start, -1, result);
      
      // conversion
      Counter<UnorderedPair<TreeNode,TreeNode>> convertedResult = new Counter<UnorderedPair<TreeNode,TreeNode>>();
      for (int l1 = 0; l1 < nLeaves; l1++)
      {
        final TreeNode t1 = indexer.i2o(l1);
        for (int l2 = l1 + 1; l2 < nLeaves; l2++)
          convertedResult.setCount(new UnorderedPair<TreeNode, TreeNode>(t1, indexer.i2o(l2)), result[l1][l2]);
      }
      return convertedResult;
    }

    private void _dfsTotalBL(double parentLen, int start, int current, int parent, double[][] result)
    {
      final int [] thisNbhrs = nbhrs[current];
      final double [] thisBLs = bls[current];
      if (thisNbhrs.length != thisBLs.length)
        throw new RuntimeException();
      for (int nIndex = 0; nIndex < thisNbhrs.length; nIndex++)
      {
        int nbhr = thisNbhrs[nIndex];
        if (nbhr != parent)
        {
          final double newLen = parentLen + thisBLs[nIndex];
          if (nbhr < nLeaves)
            // a leaf!
            result[start][nbhr] = newLen;
          else
            // recurse!
            _dfsTotalBL(newLen, start, nbhr, current, result);
        }
      }
    }
  }
  
  /**
   * Remove nodes that have exactly two neighbors, changing branch lengths appropriately 
   * @param t
   */
  public static void simplify(UnrootedTree t)
  {
    simplify(t, t.leaves().get(0), null);
  }
  
  /**
   * See comments in unrooted
   * 
   * @param t
   * @param parent
   * @param current
   * @param accumulatedBranchLength
   * @return
   */
  private static void simplify(UnrootedTree t, TreeNode goodNode, TreeNode parent)
  {
    outerLoop : for (UnorderedPair<TreeNode,TreeNode> neighborEdge : Lists.newArrayList(t.getTopology().edgesOf(goodNode)))
    {
      TreeNode neighborNode = GraphUtils.pickOther(neighborEdge, goodNode);
      if (neighborNode == parent)
        continue outerLoop;
      
      // follow down until we get to a good node,
      TreeNode previous = goodNode;
      TreeNode current = neighborNode;
      double totalBranchLength = 0.0;
      List<Pair<TreeNode,TreeNode>> path = Lists.newArrayList();
      do
      {
        totalBranchLength += t.getBranchLength(previous, current);
        path.add(Pair.of(previous, current));
        TreeNode next = findNext(t.getTopology().edgesOf(current), current, previous);
        previous = current;
        current = next;
      } while (previous != null && t.getTopology().edgesOf(previous).size() == 2);
      
      // if there were more than one hop to a good node, do simplifications
      TreeNode last = BriefLists.last(path).getRight();
      if (path.size() > 1)
      {
        for (int i = 0; i < path.size(); i++)
        {
          Pair<TreeNode,TreeNode> currentEdge = path.get(i);
          t.removeEdge(currentEdge.getLeft(), currentEdge.getRight());
          if (i != 0)
            t.getTopology().removeVertex(currentEdge.getLeft());
        }
        t.addEdge(goodNode, last, totalBranchLength);
      }
      
      simplify(t, last, goodNode);
    }

  }
  private static TreeNode findNext(
      Set<UnorderedPair<TreeNode, TreeNode>> edges, TreeNode center, TreeNode previous)
  {
    for (UnorderedPair<TreeNode, TreeNode> edge : edges)
    {
      TreeNode node = GraphUtils.pickOther(edge, center);
      if (node != previous)
        return node;
    }
    return null;
  }
  

}
