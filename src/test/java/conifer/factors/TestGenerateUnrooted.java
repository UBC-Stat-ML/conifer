package conifer.factors;

import java.util.List;
import java.util.Random;
import java.util.Set;

import org.junit.Assert;
import org.junit.Test;

import blang.core.RealConstant;
import blang.core.RealDistribution;
import blang.distributions.Exponential;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;
import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.UnrootedTreeUtils;
import conifer.UnrootedTreeUtils.Clade;



public class TestGenerateUnrooted
{
  @Test
  public void test()
  {
    Counter<Set<UnorderedPair<Clade, Clade>>> biparts = new Counter<Set<UnorderedPair<Clade,Clade>>>();
    Indexer<TreeNode> tipIndexer = null;
    RealDistribution exp = Exponential.distribution(new RealConstant(1.0)); 
    Random rand = new Random(1);
    List<TreeNode> leaves = TopologyUtils.syntheticTaxaList(6);
    
    final int nIters = 100000;
    for (int i = 0; i < nIters; i++)
    {
      UnrootedTree tree = NonClockTreePriorUtils.sample(rand, exp, leaves);
      
      if (tipIndexer == null)
        tipIndexer = UnrootedTreeUtils.tipIndexer(tree);
      biparts.incrementCount(UnrootedTreeUtils.bipartitions(tree, tipIndexer), 1.0);
      
    }
    biparts.normalize();
    
    double truth = 1.0/biparts.size();
    // TODO: clean with a stat test
    for (Set<UnorderedPair<Clade, Clade>> key : biparts)
      Assert.assertTrue(Math.abs(truth - biparts.getCount(key))/truth < 0.1);
  }

}
