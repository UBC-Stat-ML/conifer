package conifer;


import java.util.List;
import java.util.Random;

import org.jgrapht.Graphs;

import com.google.common.collect.Lists;

import conifer.factors.NonClockTreePrior;

import bayonet.distributions.Exponential;
import blang.variables.RealVariable;
import briefj.BriefCollections;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;



public class UnrootedTreeUtils
{
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
  
  public static void main(String [] args)
  {
    Random rand = new Random(1);
    for (int i = 0; i < 3; i++)
    {
      List<TreeNode> leaves = Lists.newArrayList();
      leaves.add(TreeNode.withLabel("A"));
      leaves.add(TreeNode.withLabel("B"));
      leaves.add(TreeNode.withLabel("C"));
      leaves.add(TreeNode.withLabel("D"));
      UnrootedTree randomTree = NonClockTreePrior.generate(rand, Exponential.on(RealVariable.real()), leaves);
      System.out.println(toNewick(randomTree));
    }
  }
}
