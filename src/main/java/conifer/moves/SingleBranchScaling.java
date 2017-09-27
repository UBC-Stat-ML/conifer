package conifer.moves;

import java.util.List;
import java.util.Random;

import blang.mcmc.internals.Callback;
import blang.mcmc.MHSampler;
import briefj.collections.UnorderedPair;

import com.google.common.collect.Lists;

import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.RandomUtils.DiscreteUniform;



public class SingleBranchScaling extends MHSampler<UnrootedTree>
{
  private static final double lambda = 2.0 * Math.log(2.0);

  @Override
  public void propose(Random rand, Callback callback)
  {
    List<UnorderedPair<TreeNode, TreeNode>> allEdges = Lists.newArrayList(variable.getTopology().edgeSet());
    final UnorderedPair<TreeNode, TreeNode> edge = DiscreteUniform.sample(allEdges, rand);
    final double oldValue = variable.getBranchLength(edge);
    double u = rand.nextDouble();
    final double m = Math.exp(lambda * (u - 0.5));
    final double newValue = m * oldValue;
    callback.setProposalLogRatio(Math.log(m));
    variable.updateBranchLength(edge, newValue);
    if (!callback.sampleAcceptance())
      variable.updateBranchLength(edge, oldValue);
  }
}
