package conifer.moves

import blang.mcmc.internals.Callback
import blang.mcmc.MHSampler
import blang.mcmc.SampledVariable
import static com.google.common.collect.Lists.*
import bayonet.distributions.Random
import conifer.UnrootedTree
import static conifer.Utils.*
import static java.lang.Math.*

class SingleBranchScaling extends MHSampler {
  @SampledVariable UnrootedTree variable
  
  override propose(Random rand, Callback callback) {
    val allEdges = newArrayList(variable.topology.edgeSet)
    val edge = sample(allEdges, rand)
    val oldValue = variable.getBranchLength(edge)
    val m = exp(2.0 * log(2.0) * (rand.nextDouble - 0.5))
    callback.setProposalLogRatio(log(m))
    variable.updateBranchLength(edge, m * oldValue)
    if (!callback.sampleAcceptance) 
      variable.updateBranchLength(edge, oldValue)
  }
}
