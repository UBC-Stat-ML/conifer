package conifer.factors;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.ejml.simple.SimpleMatrix;

import tutorialj.Tutorial;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.ctmc.CTMC;
import conifer.ctmc.JukeCantorRateMatrix;
import conifer.ctmc.RateMatrix;
import conifer.io.PhylogeneticObservationFactory;
import bayonet.distributions.Exponential;
import bayonet.marginal.UnaryFactor;
import bayonet.marginal.algo.EdgeSorter;
import bayonet.marginal.algo.SumProduct;
import bayonet.marginal.discrete.DiscreteFactorGraph;
import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;
import blang.factors.Factor;
import blang.variables.RealVariable;
import briefj.BriefCollections;
import briefj.BriefIO;



public class TreeLikelihood<R extends RateMatrix> implements Factor
{
  @FactorArgument
  public final UnrootedTree tree;
  
  @FactorComponent
  public final R rateMatrix;
  
  private final Map<TreeNode, UnaryFactor<TreeNode>> observationFactors;
  private final int nSites;
  
  
  
  private TreeLikelihood(UnrootedTree tree, R rateMatrix,
      Map<TreeNode, UnaryFactor<TreeNode>> observationFactors, int nSites)
  {
    this.tree = tree;
    this.rateMatrix = rateMatrix;
    this.observationFactors = observationFactors;
    this.nSites = nSites;
  }

  public static TreeLikelihood<JukeCantorRateMatrix> fromObservations(File inputFile)
  {
    // TODO: make this detect the type of observation
    PhylogeneticObservationFactory factory = PhylogeneticObservationFactory.nucleotidesFactory();
    Map<TreeNode, UnaryFactor<TreeNode>> leafLikelihoods = Maps.newHashMap();
    int nSites = -1;
    for (String line : BriefIO.readLines(inputFile))
      if (!line.matches("^\\s*$"))
      {
        String [] split = line.split("\\s+");
        if (split.length != 2)
          throw new RuntimeException();
        TreeNode leafNode = TreeNode.withLabel(split[0]);
        double [][] indicators = factory.site2CharacterIndicators(split[1]);
        if (nSites == -1)
          nSites = indicators.length;
        if (nSites != indicators.length)
          throw new RuntimeException();
        UnaryFactor<TreeNode> factor = DiscreteFactorGraph.createUnary(leafNode, new SimpleMatrix(indicators));
        leafLikelihoods.put(leafNode, factor);
      }
    // pick an arbitrary tree
    Random rand = new Random(1);
    List<TreeNode> leaves = Lists.newArrayList(leafLikelihoods.keySet());
    UnrootedTree randomTree = NonClockTreePrior.generate(rand, Exponential.on(RealVariable.real()), leaves);
    return new TreeLikelihood<JukeCantorRateMatrix>(randomTree, new JukeCantorRateMatrix(factory.nSymbols()), leafLikelihoods, nSites);
  }

  @Override
  public double logDensity()
  {
    SumProduct<TreeNode> sumProduct = new SumProduct<TreeNode>(buildFactorGraph());
    return sumProduct.logNormalization();
  }
  
//  TODO:
  
  /**
   * 
   */
  @Tutorial(showSource = false, showLink = true)
  public DiscreteFactorGraph<TreeNode> buildFactorGraph()
  {
    /* startRem throw new RuntimeException(); */
    
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
    
    /* endRem */
  }
  
  /* startRem  */

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
  
  /* endRem */
}
