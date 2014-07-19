package conifer.factors;

import java.util.Collection;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;

import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.UnrootedTreeUtils;
import conifer.UnrootedTreeUtils.Clade;

import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;



public class TestGenerateUnrooted
{
  @Test
  public void test()
  {
    Counter<Set<UnorderedPair<Clade, Clade>>> biparts = new Counter<Set<UnorderedPair<Clade,Clade>>>();
    Counter<RTree> rootedTopos = new Counter<RTree>();
    Indexer<TreeNode> tipIndexer = null;
    Exponential<RateParameterization> exp = Exponential.newExponential();
    Random rand = new Random(1);
    List<TreeNode> leaves = TopologyUtils.syntheticTaxaList(4);
    
    final int nIters = 100000;
    for (int i = 0; i < nIters; i++)
    {
      UnrootedTree tree = NonClockTreePrior.generate(rand, exp, leaves);
      
      if (tipIndexer == null)
        tipIndexer = UnrootedTreeUtils.tipIndexer(tree);
      biparts.incrementCount(UnrootedTreeUtils.bipartitions(tree, tipIndexer), 1.0);
      
      RTree topo = generate(rand, leaves);
      rootedTopos.incrementCount(topo, 1.0);
    }
    biparts.normalize();
    for (Set<UnorderedPair<Clade, Clade>> key : biparts)
      System.out.println("" + biparts.getCount(key) + "\t" + key);
    
    for (RTree topo : rootedTopos)
      System.out.println("" + rootedTopos.getCount(topo) + "\t" + topo);
  }
  
  public static class RTree
  {
    public Set<RTree> children = Sets.newLinkedHashSet();
    public TreeNode label = null;
    @Override
    public String toString()
    {
      if (label != null && !children.isEmpty())
        throw new RuntimeException();
      if (label != null)
        return label.toString();
      else
      {
        StringBuilder result = new StringBuilder();
        result.append("[");
        for (RTree child : children)
          result.append(child.toString() + " ");
        result.append("]");
        return result.toString();
      }
    }
    @Override
    public int hashCode()
    {
      final int prime = 31;
      int result = 1;
      result = prime * result + ((children == null) ? 0 : children.hashCode());
      result = prime * result + ((label == null) ? 0 : label.hashCode());
      return result;
    }
    @Override
    public boolean equals(Object obj)
    {
      if (this == obj)
        return true;
      if (obj == null)
        return false;
      if (getClass() != obj.getClass())
        return false;
      RTree other = (RTree) obj;
      if (children == null)
      {
        if (other.children != null)
          return false;
      } else if (!children.equals(other.children))
        return false;
      if (label == null)
      {
        if (other.label != null)
          return false;
      } else if (!label.equals(other.label))
        return false;
      return true;
    }
    
  }
  
  public static RTree generate(
      Random random, 
      Collection<TreeNode> leaves)
  {
    List<RTree> roots = Lists.newArrayList();
    for (TreeNode leaf : leaves)
    {
      RTree l = new RTree();
      l.label = leaf;
      roots.add(l);
    }
    while (roots.size() > 1)
    {
      
        Pair<RTree,RTree> pair = NonClockTreePrior.popRandomPair(random, roots);
        RTree parent = new RTree();
        parent.children.add(pair.getLeft());
        parent.children.add(pair.getRight());
        roots.add(parent);
    }
    return roots.get(0);
  }
}
