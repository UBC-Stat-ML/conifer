package conifer.factors;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Queue;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.math.SamplingUtils;
import blang.core.RealDistribution;
import briefj.collections.UnorderedPair;

import com.google.common.collect.Lists;

import conifer.TreeNode;
import conifer.UnrootedTree;



public class NonClockTreePriorUtils<P>
{
  /**
   * Generate a tree with a topology uniformly distributed
   * over bifurcating unrooted tree, using:
   * 
   *   The generation of random, binary unordered trees
   *   George W. Furnas
   *     see 2.1.2 (p.204)
   */
  public static UnrootedTree sample(
      Random random, 
      RealDistribution branchDistribution,
      Collection<TreeNode> leaves)
  {
    UnrootedTree result = new UnrootedTree();
    
    List<TreeNode> shuffled = Lists.newArrayList(leaves);
    Collections.shuffle(shuffled, random);
    Queue<TreeNode> queue = Lists.newLinkedList(shuffled);
    
    if (queue.isEmpty())
      return result;
    
    TreeNode leaf1 = queue.poll();
    result.addNode(leaf1);
    
    if (queue.isEmpty())
      return result;
    
    TreeNode leaf2 = queue.poll();
    result.addNode(leaf2);
    result.addEdge(leaf1, leaf2, Double.NaN);
    
    while (!queue.isEmpty())
    {
      // pick a random edge
      UnorderedPair<TreeNode, TreeNode> edge = SamplingUtils.uniformFromCollection(random, result.getTopology().edgeSet());
      TreeNode internal = TreeNode.nextUnlabelled();
      TreeNode newLeaf = queue.poll();
      result.removeEdge(edge);
      result.addNode(newLeaf);
      result.addNode(internal);
      result.addEdge(newLeaf, internal, Double.NaN);
      result.addEdge(internal, edge.getFirst(), Double.NaN);
      result.addEdge(internal, edge.getSecond(), Double.NaN);
    }
    
    for (UnorderedPair<TreeNode, TreeNode> edge : result.getTopology().edgeSet())
      result.updateBranchLength(edge, branchDistribution.sample(random));
    
    return result;
  }
  
  public static <T> Pair<T,T> popRandomPair(Random rand, List<T> items)
  {
    if (items.size() < 2) 
      throw new RuntimeException();
    List<Integer> indices = sampleWithoutReplacement(rand, items.size(), 2);
    Collections.sort(indices);
    Pair<T,T> result = Pair.of(items.get(indices.get(0)), items.get(indices.get(1)));
    items.remove((int)indices.get(1)); // remove second before first to avoid shifts
    items.remove((int)indices.get(0));
    return result;
  }
  
  /**
   * Returns n samples from {0, ..., s-1}, without replacement
   * 
   * TODO: could be more efficient
   * TODO: move to Bayonet
   * 
   * @param n
   * @param s
   * @return
   */
  public static List<Integer> sampleWithoutReplacement(Random rand, int s, int n)
  {
    if (n > s || s < 0 || n < 0)
      throw new RuntimeException();
    List<Integer> list = new ArrayList<Integer>(s),
                  result=new ArrayList<Integer>(n);
    for (int i = 0; i < s; i++)
      list.add(i);
    Collections.shuffle(list,rand);
    for (int i = 0; i < n; i++)
      result.add(list.get(i));
    return result;
  }

public P branchDistributionParameters;

}
