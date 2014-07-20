package conifer.moves;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.jgrapht.Graphs;

import bayonet.distributions.DiscreteUniform;
import bayonet.distributions.Multinomial;
import bayonet.graphs.GraphUtils;
import bayonet.marginal.UnaryFactor;
import bayonet.marginal.algo.SumProduct;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.NodeMove;
import blang.mcmc.SampledVariable;
import briefj.collections.UnorderedPair;

import com.google.common.collect.Lists;

import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.EvolutionaryModel;
import conifer.models.EvolutionaryModelUtils;
import conifer.models.LikelihoodComputationContext;


/**
 * A subtree prune regraft move.
 * 
 * See comments in execute for detailed instruction.
 * 
 * Assumes arity of 3 for internal nodes. 
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class SPRMove extends NodeMove
{
  @SampledVariable UnrootedTree tree;
  
  @ConnectedFactor NonClockTreePrior<?> treePrior;
  
  @ConnectedFactor UnrootedTreeLikelihood<?> treeLikelihood;
  
  @SuppressWarnings("unchecked")
  @Override
  public void execute(Random rand)
  {
    // nothing interesting to do if tree is only a single branch
    List<TreeNode> internalNodes = GraphUtils.internalNodes(tree.getTopology());
    if (internalNodes.isEmpty())
      return;
    
    // pick an internal node at random
    TreeNode removedRoot = DiscreteUniform.sample(internalNodes, rand);
    
    // one neighbor will from the edge to disconnect, and another one, the new root
    List<UnorderedPair<TreeNode,TreeNode>> neighbors = Lists.newArrayList(tree.getTopology().edgesOf(removedRoot));
    if (neighbors.size() != 3)
      throw new RuntimeException("Currently supporting only internal arities of 3.");
    Collections.shuffle(neighbors, rand);
    TreeNode 
      newRoot  = GraphUtils.pickOther(neighbors.get(0), removedRoot),
      detached = GraphUtils.pickOther(neighbors.get(1), removedRoot),
      third    = GraphUtils.pickOther(neighbors.get(2), removedRoot);
    
    // disconnect the tree
    UnrootedTree prunedSubtree = tree.prune(removedRoot, detached);
    
    final double 
      remRootThirdLen  =  tree.getBranchLength(removedRoot, third),
      newRootRemRootLen = tree.getBranchLength(removedRoot, newRoot);
      
    if (remRootThirdLen == 0.0 || newRootRemRootLen == 0.0)
      throw new RuntimeException("SPR does not support zero branch lengths.");
    
    // compute the factor graphs (one per potentical category) for the subtree, and run the sum product on these
    // TODO: consider more than one stem lengths
    EvolutionaryModel evolutionaryModel = treeLikelihood.evolutionaryModel;
    List<UnaryFactor<TreeNode>> prunedSubtreeMarginals = EvolutionaryModelUtils.getRootMarginalsFromFactorGraphs(EvolutionaryModelUtils.buildFactorGraphs(evolutionaryModel, prunedSubtree, removedRoot, treeLikelihood.observations), removedRoot);
    
    List<TreeNode> attachmentPoints = tree.addAuxiliaryInternalNodes(rand, removedRoot);//smallRatio, largerRatio, newRoot);
    
    // run the sum product on the main tree
    List<SumProduct<TreeNode>> mainTreeSumProducts = EvolutionaryModelUtils.getSumProductsFromFactorGraphs(EvolutionaryModelUtils.buildFactorGraphs(evolutionaryModel, tree, newRoot, treeLikelihood.observations, false), newRoot);
    
    // consider all re-attachments
    double [] samplingArray = new double[attachmentPoints.size()];
    for (int i = 0; i < samplingArray.length; i++)
    {
      TreeNode attachmentPoint = attachmentPoints.get(i);
      List<UnaryFactor<TreeNode>> currentFullUnaries = Lists.newArrayList();
      
      for (int f = 0; f < mainTreeSumProducts.size(); f++)
      {
        SumProduct<TreeNode> currentMainTreeSumProduct = mainTreeSumProducts.get(f);
        UnaryFactor<TreeNode> currentMainTreeMarginal = currentMainTreeSumProduct.computeMarginal(attachmentPoint);
        UnaryFactor<TreeNode> prunedSubtreeMarginal = prunedSubtreeMarginals.get(f);
        currentFullUnaries.add(currentMainTreeSumProduct.getFactorGraph().factorOperations().pointwiseProduct(Arrays.asList(prunedSubtreeMarginal, currentMainTreeMarginal)));
      }
      
      Pair<Double,Double> neighborBLs = neighborBranchLengths(attachmentPoint, tree);
      final double 
        b1 = neighborBLs.getLeft(), 
        b2 = neighborBLs.getRight();

      double priorFactor = 
        + treePrior.branchLengthLogDensity(b1) + treePrior.branchLengthLogDensity(b2) 
        - treePrior.branchLengthLogDensity(b1 + b2);
      
      LikelihoodComputationContext context = new LikelihoodComputationContext(currentFullUnaries);
      final double likelihood = evolutionaryModel.computeLogLikelihood(context);
      
      samplingArray[i] = priorFactor + likelihood;
    }
    
    // sample a re-attachment point
    Multinomial.expNormalize(samplingArray);
    int sampledIndex = Multinomial.sampleMultinomial(rand, samplingArray);
    
    double 
      proposedBranches = sum(neighborBranchLengths(attachmentPoints.get(sampledIndex), tree)),
      oriBranches      = sum(neighborBranchLengths(removedRoot, tree));
    double ratio = (proposedBranches) / (oriBranches);
    
    TreeNode sampledAttachment = rand.nextDouble() < ratio ? attachmentPoints.get(sampledIndex) : removedRoot;
    
    // do the re-attachment
    tree.regraft(prunedSubtree, removedRoot, sampledAttachment);
    tree.simplify();
  }

  private static double sum(Pair<Double,Double> pair)
  {
    return pair.getLeft() + pair.getRight();
  }
  
  private static Pair<Double,Double> neighborBranchLengths(TreeNode treeNode, UnrootedTree tree)
  {
    List<TreeNode> attachmentNeighbors = Graphs.neighborListOf(tree.getTopology(), treeNode);
    if (attachmentNeighbors.size() != 2)
      throw new RuntimeException();
    double
      b1 = tree.getBranchLength(attachmentNeighbors.get(0), treeNode), 
      b2 = tree.getBranchLength(attachmentNeighbors.get(1), treeNode);
    return Pair.of(b1, b2);
  }
}
