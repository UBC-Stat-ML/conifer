package conifer.moves;

import java.util.List;
import java.util.Random;

import org.jgrapht.Graphs;

import blang.mcmc.Callback;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.MHSampler;
import blang.mcmc.SampledVariable;
import briefj.collections.UnorderedPair;

import com.google.common.collect.Lists;

import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;



public class SingleNNI extends MHSampler<UnrootedTree>
{
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
      variable.interchange(moved1, referenceEdge.getFirst(), moved2, referenceEdge.getSecond());

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
