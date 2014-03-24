package conifer.models;


import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;



import org.ejml.simple.SimpleMatrix;
import org.jgrapht.UndirectedGraph;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.ctmc.CTMC;
import conifer.ctmc.RateMatrixToEmissionModel;
import conifer.factors.NonClockTreePrior;
import conifer.io.PhylogeneticObservationFactory;
import conifer.io.TreeObservations;
import bayonet.distributions.Exponential;
import bayonet.distributions.Multinomial;
import bayonet.marginal.FactorGraph;
import bayonet.marginal.UnaryFactor;
import bayonet.marginal.algo.ExactSampler;
import bayonet.marginal.discrete.DiscreteFactorGraph;
import bayonet.math.NumericalUtils;
import blang.annotations.FactorComponent;
import briefj.collections.UnorderedPair;


/**
 * A substitution-only phylogenetic likelihood model where each 
 * site can have an independent latent variable determining which
 * of a finite collection of rate matrices to use.
 * 
 * This includes multiple rate models such as Yang, J. Mol Evo 94,
 * but can also be more general.
 * 
 * This also supports a transition matrix at each leaf to model 
 * for example measurement error models, or more general setups where
 * the characters evolved along the topology come from a different
 * space thant the space of observations (e.g. covarion models).
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class MultiCategorySubstitutionModel implements EvolutionaryModel
{
  @FactorComponent
  public final RateMatrixMixture rateMatrixMixture;
  
  public final int nSites;
  
  public MultiCategorySubstitutionModel(RateMatrixMixture rateMatrixMixture, int nSites)
  {
    this.rateMatrixMixture = rateMatrixMixture;
    this.nSites = nSites;
  }

  @Override
  public int nFactorGraphs()
  {
    return RateMatrixMixtureUtils.nCategories(rateMatrixMixture);
  }

  @Override
  public FactorGraph<TreeNode> newFactorGraph(
      UndirectedGraph<TreeNode, UnorderedPair<TreeNode, TreeNode>> undirectedGraph)
  {
    return new DiscreteFactorGraph<TreeNode>(undirectedGraph);
  }

  @Override
  public void buildObservation(TreeNode leaf, LikelihoodFactoryContext context)
  {
    final int categoryIndex = context.getFactorGraphIndex();
    double [][] observation = (double[][]) context.getObservations().get(leaf);
    UnaryFactor<TreeNode> observationUnary = DiscreteFactorGraph.createUnary(observation);
    RateMatrixToEmissionModel emissionModel = rateMatrixMixture.getRateMatrix(categoryIndex).getEmissionModel();
    if (emissionModel != null)
    {
      double [][] latent2Observation = emissionModel.getMatrixStatesToObservationProbabilities();
      observationUnary = context.getDiscreteFactorGraph().marginalize(new SimpleMatrix(latent2Observation), observationUnary);
    }
    context.getDiscreteFactorGraph().unaryTimesEqual(leaf, observationUnary);
  }

  @Override
  public void buildTransition(TreeNode parent, TreeNode children,
      LikelihoodFactoryContext context)
  {
    final int categoryIndex = context.getFactorGraphIndex();
    CTMC ctmc = getCTMC(context.getCache(), categoryIndex);
    double branchLength = context.getBranchLength(parent, children);
    double [][] transitionMatrix = ctmc.marginalTransitionProbability(branchLength);
    context.getDiscreteFactorGraph().setBinary(parent, children, new SimpleMatrix(transitionMatrix));
  }

  private final Object CTMC_KEY = new Object();
  @SuppressWarnings({ "rawtypes", "unchecked" })
  private CTMC getCTMC(Map cache, int categoryIndex)
  {
    if (cache.containsKey(CTMC_KEY))
      return (CTMC) cache.get(CTMC_KEY);
    CTMC result = rateMatrixMixture.getRateMatrix(categoryIndex).getProcess(); //new CTMC(rateMatrixMixture.getRateMatrix(categoryIndex));
    cache.put(CTMC_KEY, result);
    return result;
  }

  @Override
  public void buildInitialDistribution(TreeNode node,
      LikelihoodFactoryContext context)
  {
    final int categoryIndex = context.getFactorGraphIndex();
    CTMC ctmc = getCTMC(context.getCache(), categoryIndex);
    double [] stationary = ctmc.stationaryDistribution();
    double [][] allStatios = new double[nSites][stationary.length];
    for (int siteIndex = 0; siteIndex < nSites; siteIndex++)
      allStatios[siteIndex] = stationary;
    context.getDiscreteFactorGraph().unaryTimesEqual(node, new SimpleMatrix(allStatios));
  }

  @Override
  public double computeLogLikelihood(LikelihoodComputationContext context)
  {
    // Note: not particularly efficient: lots of logs and array accessed in bad ways,
    // but this occurs only at one point of the tree, so should not be a huge bottleneck in large trees
    // (but could be when the number of sites is large; then could rewrite this with scalings)
    List<UnaryFactor<TreeNode>> rootMarginals = context.getRootMarginals();
    final int nCat = rootMarginals.size();
    List<Double> categoryPriorLogPrs = rateMatrixMixture.getLogPriorProbabilities(); 
    double[][] categoryAndSiteSpecificLikelihoods = new double[nCat][];
    for (int c = 0; c < rootMarginals.size(); c++)
      categoryAndSiteSpecificLikelihoods[c] = DiscreteFactorGraph.siteLogNormalizations(rootMarginals.get(c));
    final int nSites = categoryAndSiteSpecificLikelihoods[0].length;
    final double [] workArray = new double[nCat];
    double sum = 0.0;
    for (int s = 0; s < nSites; s++)
    {
      for (int c = 0; c < nCat; c++)
        workArray[c] = categoryAndSiteSpecificLikelihoods[c][s] + categoryPriorLogPrs.get(c);
      sum += NumericalUtils.logAdd(workArray);
    }
    return sum;
  }

  @SuppressWarnings({ "unchecked", "rawtypes" })
  @Override
  public void generateObservationsInPlace(
      Random rand, 
      TreeObservations destination,
      UnrootedTree tree, TreeNode root)
  {
    // sample full paths for all categories (a bit wasteful, but not more costly than doing posterior inference)
    List<FactorGraph<TreeNode>> factorGraphs = EvolutionaryModelUtils.buildFactorGraphs(this, tree, root, destination);
    List<Map<TreeNode, double[][]>> allSamples = Lists.newArrayList();
    int nSites = -1;
    for (FactorGraph<TreeNode> factorGraph : factorGraphs)
    {
      ExactSampler<TreeNode> sampler = ExactSampler.priorSampler(factorGraph, ((DiscreteFactorGraph)factorGraph).getSampler());
      Map<TreeNode, UnaryFactor<TreeNode>> samples = sampler.sample(rand, root);
      Map<TreeNode, double[][]> transformedSamples = Maps.newHashMap();
      for (TreeNode leaf : TopologyUtils.leaves(tree.getTopology()))
      {
        double [][] current = DiscreteFactorGraph.getNormalizedCopy(samples.get(leaf));
        nSites = current.length;
        transformedSamples.put(leaf, current);
      }
      allSamples.add(transformedSamples);
    }
    // sample category indicators
    final int nCat = factorGraphs.size();
    int [] indicators = new int[nCat];
    double [] prs = RateMatrixMixtureUtils.priorProbabilities(rateMatrixMixture);
    for (int site = 0; site < nSites; site++)
      indicators[site] = Multinomial.sampleMultinomial(rand, prs);
    // compile the result
    for (TreeNode leaf : TopologyUtils.leaves(tree.getTopology()))
    {
      double [][][] allSamplesTransformed = new double[nCat][][]; // cat -> site -> state
      for (int cat = 0; cat < nCat; cat++)
        allSamplesTransformed[cat] = allSamples.get(cat).get(leaf);
      
      double [][] result = new double[nSites][]; // site -> state
      for (int site = 0; site < nSites; site++)
        result[site] = allSamplesTransformed[indicators[site]][site];
      
      destination.set(leaf, result);
    }
  }

  public static void loadObservations(
      TreeObservations observations, 
      LinkedHashMap<TreeNode, CharSequence> rawStrings, 
      PhylogeneticObservationFactory observationFactory)
  {
    for (TreeNode treeNode : rawStrings.keySet())
    {
      String rawString = rawStrings.get(treeNode).toString();
      double [][] indicators = observationFactory.site2CharacterIndicators(rawString);
      observations.set(treeNode, indicators);
    }
  }
  
  /**
   * 
   * @param leaves
   * @return A default tree obtained by sampling a nonclock tree with an exponential rate one
   *         and a fixed seed.
   */
  public static UnrootedTree defaultTree(List<TreeNode> leaves)
  {
    Random rand = new Random(1);
    return NonClockTreePrior.generate(rand, Exponential.newExponential(), leaves);
  }

}