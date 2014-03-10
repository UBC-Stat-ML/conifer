package conifer;

import java.util.List;
import java.util.Random;

import org.jgrapht.Graphs;

import com.google.common.collect.Lists;


import conifer.factors.NonClockTreePrior;

import bayonet.distributions.Exponential;
import blang.variables.RealVariable;
import briefj.BriefCollections;



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
