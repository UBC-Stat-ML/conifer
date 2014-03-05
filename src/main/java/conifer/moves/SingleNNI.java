package conifer.moves;

import java.util.List;
import java.util.Random;

import org.jgrapht.Graphs;

import bayonet.distributions.DiscreteUniform;
import blang.factors.Factor;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.MHProposalDistribution;
import blang.mcmc.SampledVariable;

import com.google.common.collect.Lists;

import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
import fig.basic.UnorderedPair;



public class SingleNNI implements MHProposalDistribution
{
  // TODO: use constructor signature instead
  // rationale: removes the unused warning issue
  // allows for different behavior in the same compilation unit
  // note: still use annotations in the constructor signature?
  // con: makes code more verbose
  
  @SampledVariable UnrootedTree tree;
  
  @ConnectedFactor List<Factor> factors;

  @Override
  public Proposal propose(Random rand)
  {
    List<UnorderedPair<TreeNode, TreeNode>> nonTerminalEdges = TopologyUtils.nonTerminalEdges(tree.getTopology());
    UnorderedPair<TreeNode, TreeNode> referenceEdge = DiscreteUniform.sample(nonTerminalEdges, rand);
    TreeNode moved1 = sampleMovedEndpoint(referenceEdge, referenceEdge.getFirst(), rand);
    TreeNode moved2 = sampleMovedEndpoint(referenceEdge, referenceEdge.getSecond(), rand);
    tree.interchange(moved1, referenceEdge.getFirst(), moved2, referenceEdge.getSecond());
    return new ProposalRealization(moved2, referenceEdge.getFirst(), moved1, referenceEdge.getSecond());
  }
  
  private TreeNode sampleMovedEndpoint(
      UnorderedPair<TreeNode, TreeNode> referenceEdge, 
      TreeNode fixedEndPoint,
      Random rand)
  {
    List<TreeNode> choices = Lists.newArrayList();
    for (TreeNode neibhbor : Graphs.neighborListOf(tree.getTopology(), fixedEndPoint))
      if (!neibhbor.equals(referenceEdge.getFirst()) &&
          !neibhbor.equals(referenceEdge.getSecond()))
        choices.add(neibhbor);
    if (choices.size() < 2)
      throw new RuntimeException();
    return DiscreteUniform.sample(choices, rand);
  }
  
  private class ProposalRealization implements Proposal
  {
    private final TreeNode moved1, fixed1, moved2, fixed2;
    private boolean done = false;
    
    private ProposalRealization(TreeNode moved1, TreeNode fixed1,
        TreeNode moved2, TreeNode fixed2)
    {
      this.moved1 = moved1;
      this.fixed1 = fixed1;
      this.moved2 = moved2;
      this.fixed2 = fixed2;
    }
    @Override
    public double logProposalRatio()
    {
      return 0;
    }
    @Override
    public void acceptReject(boolean accept)
    {
      if (done)
        throw new RuntimeException();
      done = true;
      if (!accept)
        tree.interchange(moved1, fixed1, moved2, fixed2);
    }
  }

}
