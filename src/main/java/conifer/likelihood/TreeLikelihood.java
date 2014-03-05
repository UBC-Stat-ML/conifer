package conifer.likelihood;

import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.ejml.simple.SimpleMatrix;

import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.ctmc.CTMC;
import conifer.ctmc.RateMatrix;
import bayonet.marginal.UnaryFactor;
import bayonet.marginal.algo.EdgeSorter;
import bayonet.marginal.algo.SumProduct;
import bayonet.marginal.discrete.DiscreteFactorGraph;
import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;
import blang.factors.Factor;
import briefj.BriefCollections;



public class TreeLikelihood<R extends RateMatrix> implements Factor
{
  @FactorArgument
  private final UnrootedTree tree;
  
  @FactorComponent
  private final R rateMatrix;
  
  private final Map<TreeNode, UnaryFactor<TreeNode>> observationFactors;
  private final int nSites;
  
  
  public static <R extends RateMatrix> TreeLikelihood<R> fromObservations(File inputFile)
  {
    
  }

  @Override
  public double logDensity()
  {
    SumProduct<TreeNode> sumProduct = new SumProduct<TreeNode>(buildFactorGraph());
    return sumProduct.logNormalization();
  }
  
  public DiscreteFactorGraph<TreeNode> buildFactorGraph()
  {
    // graph
    DiscreteFactorGraph<TreeNode> result = new DiscreteFactorGraph<TreeNode>(tree.getTopology());
    
    // orient edges away from root
    TreeNode arbitraryRoot = BriefCollections.pick(tree.getTopology().vertexSet());
    List<Pair<TreeNode,TreeNode>> orientedEdges = EdgeSorter.newEdgeSorter(tree.getTopology(), arbitraryRoot).backwardMessages();
    
    // observations
    addObservationFactors(result);
    
    // build CTMC from rate matrix
    CTMC ctmc = new CTMC(rateMatrix);
    
    // initial distribution
    addInitialDistribution(result, arbitraryRoot, nSites, ctmc);
    
    // transition
    addTransitions(result, orientedEdges, ctmc);
    
    return result;
  }

  private void addTransitions(DiscreteFactorGraph<TreeNode> result, List<Pair<TreeNode, TreeNode>> orientedEdges, CTMC ctmc)
  {
    for (Pair<TreeNode,TreeNode> edge : orientedEdges)
    {
      double branchLength = tree.getBranchLength(edge.getLeft(), edge.getRight());
      double [][] transitionMatrix = ctmc.marginalTransitionProbability(branchLength);
      result.setBinary(edge.getLeft(), edge.getRight(), new SimpleMatrix(transitionMatrix));
    }
  }

  private void addInitialDistribution(DiscreteFactorGraph<TreeNode> result, TreeNode arbitraryRoot, int nSites, CTMC ctmc)
  {
    double [] stationary = ctmc.stationaryDistribution();
    double [][] allStatios = new double[nSites][stationary.length];
    for (int siteIndex = 0; siteIndex < nSites; siteIndex++)
      allStatios[siteIndex] = stationary;
    result.unariesTimesEqual(arbitraryRoot, new SimpleMatrix(allStatios));
  }

  private void addObservationFactors(DiscreteFactorGraph<TreeNode> result)
  {
    for (TreeNode observedNode : observationFactors.keySet())
      result.setUnary(observedNode, observationFactors.get(observedNode));
  }
  
  
}
