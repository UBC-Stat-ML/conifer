package conifer.moves;

import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import com.google.common.collect.Lists;

import conifer.TreeNode;
import conifer.UnrootedTree;

import bayonet.distributions.DiscreteUniform;
import bayonet.distributions.Exponential;
import blang.factors.Factor;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.MHProposalDistribution;
import blang.mcmc.SampledVariable;
import briefj.BriefLog;
import briefj.collections.UnorderedPair;



public class SingleBranchScaling implements MHProposalDistribution
{
  @SampledVariable UnrootedTree tree;
  
  @ConnectedFactor List<Factor> factors;
  
  private static final double lambda = 2.0 * Math.log(2.0);

  @Override
  public Proposal propose(Random rand)
  {
//    BriefLog.warnOnce("disabled SBS");
//    if (true)
//      return null;
    
    List<UnorderedPair<TreeNode, TreeNode>> allEdges = Lists.newArrayList(tree.getTopology().edgeSet());
    final UnorderedPair<TreeNode, TreeNode> edge = DiscreteUniform.sample(allEdges, rand);
    final double oldValue = tree.getBranchLength(edge);
    double u = rand.nextDouble();
    final double m = Math.exp(lambda * (u - 0.5));
    final double newValue = m * oldValue;
    
    tree.updateBranchLength(edge, newValue);
    return new Proposal() {
      
      @Override
      public double logProposalRatio()
      {
        return Math.log(m);
      }
      
      @Override
      public void acceptReject(boolean accept)
      {
        if (!accept)
          tree.updateBranchLength(edge, oldValue);
      }
    };
  }
  
  public static void main(String [] args)
  {
    Random rand = new Random(1);
    for (int j = 0; j < 100; j++)
    {
      
      double current = Exponential.generate(rand, 1.0);
      SummaryStatistics stat = new SummaryStatistics();
      for (int i = 0; i <10000; i++)
      {
        double oldValue = current;
        final double u = rand.nextDouble();
        final double m = Math.exp(lambda * (u - 0.5));
        final double newValue = m * oldValue;
        
        final double logRatio =  Math.log(m) 
          + Exponential.logDensity(newValue, 1.0) 
          - Exponential.logDensity(oldValue, 1.0);
        
        final boolean accept = rand.nextDouble() < Math.exp(logRatio);
          
        if (accept)
          current = newValue;
        else
          current = oldValue;
        
        stat.addValue(current);
  
      }
      System.out.println(stat.getMean());
    }
  }
}
