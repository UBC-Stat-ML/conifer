package conifer.moves;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jgrapht.Graphs;

import com.google.common.collect.Lists;

import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.UnrootedTreeUtils;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.EvolutionaryModel;
import conifer.models.EvolutionaryModelUtils;
import conifer.models.LikelihoodComputationContext;

import bayonet.distributions.DiscreteUniform;
import bayonet.distributions.Exponential;
import bayonet.distributions.Multinomial;
import bayonet.graphs.GraphUtils;
import bayonet.marginal.UnaryFactor;
import bayonet.marginal.algo.SumProduct;
import bayonet.math.NumericalUtils;
import bayonet.math.SamplingUtils;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.NodeMove;
import blang.mcmc.SampledVariable;
import briefj.BriefLog;
import briefj.collections.UnorderedPair;


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
  
  private SummaryStatistics stat = new SummaryStatistics();
  
//  @SuppressWarnings("unchecked")
//  @Override
//  public void execute(Random rand)
//  {
//    double oriLL = treeLikelihood.logDensity();
//    UnrootedTree initial = new UnrootedTree(tree);
//    
//    // nothing interesting to do if tree is only a single branch
//    List<TreeNode> internalNodes = GraphUtils.internalNodes(tree.getTopology());
//    if (internalNodes.isEmpty())
//      return;
//    
//    TreeNode removedRoot = DiscreteUniform.sample(internalNodes, rand);
//    
//    // one neighbor will from the edge to disconnect, and another one, the new root
//    List<UnorderedPair<TreeNode,TreeNode>> neighbors = Lists.newArrayList(tree.getTopology().edgesOf(removedRoot));
//    if (neighbors.size() != 3)
//      throw new RuntimeException("Currently supporting only internal arities of 3.");
//    Collections.shuffle(neighbors, rand);
//    TreeNode 
//      newRoot  = GraphUtils.pickOther(neighbors.get(0), removedRoot),
//      detached = GraphUtils.pickOther(neighbors.get(1), removedRoot),
//      third    = GraphUtils.pickOther(neighbors.get(2), removedRoot);
//    
//    final double 
//    remRootThirdLen = tree.getBranchLength(removedRoot, third),
//    newRootRemRootLen=tree.getBranchLength(removedRoot, newRoot);
//  
//    // calculate the branch ratio
//    double referenceLengthFromBot = remRootThirdLen;
//    double joinedLength = referenceLengthFromBot + newRootRemRootLen;
//    
//    double randomRatio = rand.nextDouble();
//    double b1 = 
//      randomRatio * joinedLength;
////      Exponential.generate(rand, 1.0);
//    double b2 = 
//      (1.0 - randomRatio) * joinedLength;
////      Exponential.generate(rand, 1.0);
//    tree.updateBranchLength(UnorderedPair.of(newRoot, removedRoot), b1); //
//    tree.updateBranchLength(UnorderedPair.of(removedRoot, third), b2); //
//    
//    double newLL = treeLikelihood.logDensity();
//    double logRatio = newLL - oriLL;
//    
//    if (rand.nextDouble() < Math.exp(logRatio)) // accept
//      ;
////      stat.addValue(1.0);
//    else
//    {
////      stat.addValue(0.0);
//      tree.setTo(initial);
////      System.out.println("ref, ll = " + treeLikelihood.logDensity());
//    }
////    if (stat.getN() % 10 == 0.0)
////      System.out.println(stat.getMean());
//  }
  
//  @SuppressWarnings("unchecked")
//  @Override
//  public void execute(Random rand)
//  {
//    
//    double oriLL = treeLikelihood.logDensity();
//    System.out.println("---");
//    System.out.println("ori LL = " + oriLL);
//    UnrootedTree initial = new UnrootedTree(tree);
//    
//    // nothing interesting to do if tree is only a single branch
//    List<TreeNode> internalNodes = GraphUtils.internalNodes(tree.getTopology());
//    if (internalNodes.isEmpty())
//      return;
//    
//    // pick an internal node at random
//    TreeNode removedRoot = DiscreteUniform.sample(internalNodes, rand);
//    
//    // one neighbor will from the edge to disconnect, and another one, the new root
//    List<UnorderedPair<TreeNode,TreeNode>> neighbors = Lists.newArrayList(tree.getTopology().edgesOf(removedRoot));
//    if (neighbors.size() != 3)
//      throw new RuntimeException("Currently supporting only internal arities of 3.");
//    Collections.shuffle(neighbors, rand);
//    TreeNode 
//      newRoot  = GraphUtils.pickOther(neighbors.get(0), removedRoot),
//      detached = GraphUtils.pickOther(neighbors.get(1), removedRoot),
//      third    = GraphUtils.pickOther(neighbors.get(2), removedRoot);
//    
//    // disconnect the tree
//    UnrootedTree prunedSubtree = tree.prune(removedRoot, detached);
//    
//    final double 
//      remRootThirdLen = tree.getBranchLength(removedRoot, third),
//      newRootRemRootLen=tree.getBranchLength(removedRoot, newRoot);
//    
//    // calculate the branch ratio
//    double referenceLengthFromBot = remRootThirdLen;
//    double joinedLength = referenceLengthFromBot + newRootRemRootLen;
////    double fixedRatioFromBot = referenceLengthFromBot / joinedLength;
//    
//    tree.simplify(third, removedRoot, newRoot);
//    
//    // pick a candidate edge at random
//    UnorderedPair<TreeNode, TreeNode> reattachEdge = UnorderedPair.of(newRoot, third);
//      //SamplingUtils.uniformFromCollection(rand, tree.getTopology().edgeSet());
//    double reattachBL = tree.getBranchLength(reattachEdge);
//    tree.removeEdge(reattachEdge);
//    TreeNode reattachPoint = TreeNode.nextUnlabelled();
//    tree.addNode(reattachPoint);
//    double reattachInterval = rand.nextDouble() * reattachBL;
//    tree.addEdge(reattachEdge.getFirst(), reattachPoint, reattachInterval);
//    tree.addEdge(reattachEdge.getSecond(), reattachPoint, reattachBL - reattachInterval);
//    
//    tree.regraft(prunedSubtree, removedRoot, reattachPoint);
//    
//    double newLL = treeLikelihood.logDensity();
//    double logRatio = newLL - oriLL - Math.log(joinedLength) + Math.log(reattachBL);
//    
//    System.out.println(logRatio);
//    if (rand.nextDouble() < Math.exp(logRatio)) // accept
//      System.out.println("accept!");
//    else
//    {
//      tree.setTo(initial);
//      System.out.println("ref, ll = " + treeLikelihood.logDensity());
//    }
//  }
  
  @SuppressWarnings("unchecked")
  @Override
  public void execute(Random rand)
  {
    
//    System.out.println("----");
//    System.out.println(UnrootedTreeUtils.totalTreeLength(tree));
//    System.out.println("original LL = " + treeLikelihood.logDensity());
    BriefLog.warnOnce("remove hacks in SPRMove");
    double oriLL = treeLikelihood.logDensity();
    int initIndex = -1;
    
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
      remRootThirdLen = tree.getBranchLength(removedRoot, third),
      newRootRemRootLen=tree.getBranchLength(removedRoot, newRoot);
      
    if (remRootThirdLen == 0.0 || newRootRemRootLen == 0.0)
      throw new RuntimeException("SPR does not support zero branch lengths.");
    
    // calculate the branch ratio
    double referenceLengthFromBot = remRootThirdLen;
    double joinedLength = referenceLengthFromBot + newRootRemRootLen;
    double fixedRatioFromBot = referenceLengthFromBot / joinedLength;
    
    if (fixedRatioFromBot == 1.0) // rounding problem: can cause branches of len zero
      fixedRatioFromBot = 1.0 - 1e-15; // note that double are less precise around 1.0 than around 0.0
    
    // compute the factor graphs (one per potentical category) for the subtree, and run the sum product on these
    // TODO: consider more than one stem lengths
    EvolutionaryModel evolutionaryModel = treeLikelihood.evolutionaryModel;
    List<UnaryFactor<TreeNode>> prunedSubtreeMarginals = EvolutionaryModelUtils.getRootMarginalsFromFactorGraphs(EvolutionaryModelUtils.buildFactorGraphs(evolutionaryModel, prunedSubtree, removedRoot, treeLikelihood.observations), removedRoot);
    
    // form an edge joining the two edge connected to a node of arity two caused by the disconnect
    tree.simplify(third, removedRoot, newRoot);
    
    // create intermediate nodes in the main tree
    
    
    double additionalRatio = fixedRatioFromBot; BriefLog.warnOnce("ANother thing in SPRMove to remove/fix");
      //rand.nextDouble();
    double smallRatio = Math.min(additionalRatio, fixedRatioFromBot);
    double largerRatio= Math.max(additionalRatio, fixedRatioFromBot);
    List<TreeNode> attachmentPoints = tree.addAuxiliaryInternalNodes(rand, fixedRatioFromBot, Pair.of(newRoot, third), newRoot);//smallRatio, largerRatio, newRoot);
    
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
      
      double priorFactor = 0.0;
      BriefLog.warnOnce("non exp currently not in SPR");
//      Pair<Double,Double> proposedBranches = findBranches(attachmentPoint, tree);
//      final double 
//        b1 = proposedBranches.getLeft(),
//        b2 = proposedBranches.getRight();
//      
////      System.out.println(b1 + "\t" + b2);
////      System.out.println(treePrior.branchLengthLogDensity(b1) + "\t" + treePrior.branchLengthLogDensity(b2) 
////          + "\t" + treePrior.branchLengthLogDensity(b1 + b2));
//
//      
//      double priorFactor = 
//        + treePrior.branchLengthLogDensity(b1) + treePrior.branchLengthLogDensity(b2) 
//        - treePrior.branchLengthLogDensity(b1 + b2);
      
//      System.out.println(priorFactor);
      
      LikelihoodComputationContext context = new LikelihoodComputationContext(currentFullUnaries);
      final double likelihood = evolutionaryModel.computeLogLikelihood(context);
//      System.out.println(likelihood);
      
      BriefLog.warnOnce("yet another thing to rem");
      if (NumericalUtils.isClose(likelihood, oriLL, 1e-10))
      {
//        System.out.println("*");
        initIndex = i;
      } 
      
      samplingArray[i] = priorFactor + likelihood;
    }
    
    // sample a re-attachment point
    Multinomial.expNormalize(samplingArray);
    int sampledIndex = Multinomial.sampleMultinomial(rand, samplingArray);
    
    
    BriefLog.warnOnce("yet another SPR");
    double 
      proposedBranches = findBranches2(attachmentPoints.get(sampledIndex), tree),
      oriBranches = findBranches2(attachmentPoints.get(initIndex), tree);
    double ratio = (proposedBranches) / (oriBranches);
//    System.out.println(ratio);
    if (rand.nextDouble() > ratio)   
      sampledIndex  = initIndex;
    // end bad block
    
    TreeNode sampledAttachment = attachmentPoints.get(sampledIndex);
    
    // do the re-attachment
    tree.regraft(prunedSubtree, removedRoot, sampledAttachment);
    tree.simplify();
    
//    System.out.println(treeLikelihood.logDensity());

//    System.out.println(UnrootedTreeUtils.totalTreeLength(tree));
//    System.out.println("---");
  }

  private static double findBranches2(TreeNode treeNode, UnrootedTree tree2)
  {
    double sum = 0.0;
    for (TreeNode neighbor : Graphs.neighborListOf(tree2.getTopology(), treeNode))
      sum += tree2.getBranchLength(neighbor, treeNode);
    return sum;
  }

  private static Pair<Double, Double> findBranches(TreeNode attachmentPoint,
      UnrootedTree tree)
  {
    TreeNode endPoint1 = findNeighborOfGivenDegree(attachmentPoint, tree, false);
    TreeNode attach2 = findNeighborOfGivenDegree(attachmentPoint, tree, true);
    TreeNode endPoint2 = findNeighborOfGivenDegree(attach2, tree, false);
    
    return Pair.of(tree.getBranchLength(endPoint1, attachmentPoint), tree.getBranchLength(attachmentPoint, attach2) + tree.getBranchLength(attach2, endPoint2));
  }
  
  private static TreeNode findNeighborOfGivenDegree(TreeNode node,
      UnrootedTree tree, boolean isDegree2)
  {
    TreeNode result = null;
    for (TreeNode neighbor : Graphs.neighborListOf(tree.getTopology(), node))
    {
      int currentDegree = Graphs.neighborListOf(tree.getTopology(), neighbor).size();
      boolean degreeMatch = isDegree2 ? (currentDegree == 2) : (currentDegree == 1 || currentDegree == 3);
      if (degreeMatch)
      {
        if (result != null)
          throw new RuntimeException();
        result = neighbor;
      }
    }
    if (result == null)
      throw new RuntimeException();
    return result;
  }
}
