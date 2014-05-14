package conifer.moves;

import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;

import conifer.TreeNode;
import conifer.UnrootedTree;

import bayonet.distributions.DiscreteUniform;
import blang.factors.Factor;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.MHProposalDistribution;
import blang.mcmc.SampledVariable;
import briefj.collections.UnorderedPair;



public class SingleBranchScaling implements MHProposalDistribution
{
  @SampledVariable UnrootedTree tree;
  
  @ConnectedFactor List<Factor> factors;
  
  private static final double lambda = 2.0 * Math.log(2.0);

  @Override
  public Proposal propose(Random rand)
  {
//    System.out.println("Computing single branch scaling move");
    List<UnorderedPair<TreeNode, TreeNode>> allEdges = Lists.newArrayList(tree.getTopology().edgeSet());
    final UnorderedPair<TreeNode, TreeNode> edge = DiscreteUniform.sample(allEdges, rand);
    final double oldValue = tree.getBranchLengths().get(edge);
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
  

  
  // TODO: need to merge bayonet.distribution, bayonet.mcmc.dist, fig stuff, etc
  // at same time, check all for correctness
  // TODO: add an argument to say the logDensity is normalized with respect to subset of variables
  public static double nextDouble(Random rand, double left, double right)
  {
    return left + (right-left) * rand.nextDouble();
  }

}
