package conifer.factors;

import java.util.Collection;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.junit.Assert;
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
    Indexer<TreeNode> tipIndexer = null;
    Exponential<RateParameterization> exp = Exponential.newExponential();
    Random rand = new Random(1);
    List<TreeNode> leaves = TopologyUtils.syntheticTaxaList(6);
    
    final int nIters = 100000;
    for (int i = 0; i < nIters; i++)
    {
      UnrootedTree tree = NonClockTreePrior.generate(rand, exp, leaves);
      
      if (tipIndexer == null)
        tipIndexer = UnrootedTreeUtils.tipIndexer(tree);
      biparts.incrementCount(UnrootedTreeUtils.bipartitions(tree, tipIndexer), 1.0);
      
    }
    biparts.normalize();
    
    double truth = 1.0/biparts.size();
    // TODO: clean with a stat test
    for (Set<UnorderedPair<Clade, Clade>> key : biparts)
    {
//      if (Math.abs(truth - biparts.getCount(key))/truth > 0.1)
//        System.out.println("" + truth + " vs " + biparts.getCount(key));
//      else
//        System.out.print(".");
      Assert.assertTrue(Math.abs(truth - biparts.getCount(key))/truth < 0.1);
    }
  }

}
