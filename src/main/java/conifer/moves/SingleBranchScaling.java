package conifer.moves;

import java.util.List;
import java.util.Random;

import bayonet.distributions.DiscreteUniform;
import blang.factors.Factor;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.MHProposalDistribution;
import blang.mcmc.SampledVariable;
import briefj.collections.UnorderedPair;

import com.google.common.collect.Lists;

import conifer.TreeNode;
import conifer.UnrootedTree;



public class SingleBranchScaling implements MHProposalDistribution
{
  @SampledVariable UnrootedTree tree;
  
  @ConnectedFactor List<Factor> factors;
  
  private static final double lambda = 2.0 * Math.log(2.0);

  @Override
  public Proposal propose(Random rand)
  {
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
}
