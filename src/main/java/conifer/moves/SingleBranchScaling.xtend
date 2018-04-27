package conifer.moves

import java.util.List
import blang.mcmc.internals.Callback
import blang.mcmc.MHSampler
import blang.mcmc.SampledVariable
import briefj.collections.UnorderedPair
import com.google.common.collect.Lists
import bayonet.distributions.Random
import conifer.TreeNode
import conifer.UnrootedTree
import conifer.Utils

class SingleBranchScaling extends MHSampler {
  @SampledVariable package UnrootedTree variable
  override void propose(Random rand, Callback callback) {
    var List<UnorderedPair<TreeNode, TreeNode>> allEdges = Lists.newArrayList(variable.getTopology().edgeSet())
    val UnorderedPair<TreeNode, TreeNode> edge = Utils.sample(allEdges, rand)
    val double oldValue = variable.getBranchLength(edge)
    var double u = rand.nextDouble()
    val double m = Math.exp(2.0 * Math.log(2.0) * (u - 0.5))
    val double newValue = m * oldValue
    callback.setProposalLogRatio(Math.log(m))
    variable.updateBranchLength(edge, newValue)
    if(!callback.sampleAcceptance()) variable.updateBranchLength(edge, oldValue)
  }
}
