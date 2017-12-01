package conifer.moves;

import java.util.List;

import org.jgrapht.Graphs;

import com.google.common.collect.Lists;

import bayonet.distributions.Random;
import blang.mcmc.internals.Callback;
import blang.mcmc.MHSampler;
import blang.mcmc.SampledVariable;
import briefj.collections.UnorderedPair;
import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.RandomUtils.DiscreteUniform;


public class SingleNNI extends MHSampler
{
  @SampledVariable
  UnrootedTree variable;
  
  @Override
  public void propose(Random rand, Callback callback)
  {
    List<UnorderedPair<TreeNode, TreeNode>> nonTerminalEdges = TopologyUtils.nonTerminalEdges(variable.getTopology());
    if (nonTerminalEdges.isEmpty())
      return;
    callback.setProposalLogRatio(0.0);
    UnorderedPair<TreeNode, TreeNode> referenceEdge = DiscreteUniform.sample(nonTerminalEdges, rand);
    TreeNode moved1 = sampleMovedEndpoint(referenceEdge, referenceEdge.getFirst(), rand);
    TreeNode moved2 = sampleMovedEndpoint(referenceEdge, referenceEdge.getSecond(), rand);
    variable.interchange(moved1, referenceEdge.getFirst(), moved2, referenceEdge.getSecond());
    if (!callback.sampleAcceptance())
      variable.interchange(moved2, referenceEdge.getFirst(), moved1, referenceEdge.getSecond());
  }
  
  private TreeNode sampleMovedEndpoint(
      UnorderedPair<TreeNode, TreeNode> referenceEdge, 
      TreeNode fixedEndPoint,
      Random rand)
  {
    List<TreeNode> choices = Lists.newArrayList();
    for (TreeNode neibhbor : Graphs.neighborListOf(variable.getTopology(), fixedEndPoint))
      if (!neibhbor.equals(referenceEdge.getFirst()) &&
          !neibhbor.equals(referenceEdge.getSecond()))
        choices.add(neibhbor);
    if (choices.size() < 2)
      throw new RuntimeException();
    return DiscreteUniform.sample(choices, rand);
  }

}
