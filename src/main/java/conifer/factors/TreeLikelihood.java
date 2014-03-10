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


/**
 * A factor to compute the likelihood of observed sequence data 
 * viewed as the leaves of a CTMC along a tree shaped branching 
 * process.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 * @param <R> The type of rate matrix used.
 */
public class TreeLikelihood<R extends RateMatrix> implements Factor
{
  /**
   * 
   */
  @FactorArgument
  public final UnrootedTree tree;
  
  /**
   * 
   */
  @FactorComponent
  public final R rateMatrix;
  
  private final Map<TreeNode, UnaryFactor<TreeNode>> observationFactors;
  
  /**
   * The number of sites (columns in the sequence alignment). Each is
   * considered to be independent but evolving along a shared tree.
   */
  private final int nSites;

  /**
   * Builds a factor graph, and use the elimination algorithm covered in
   * class to compute the likelihood of the observed sequences given the
   *  tree and the rate matrix.
   */
  @Override
  public double logDensity()
  {
    SumProduct<TreeNode> sumProduct = new SumProduct<TreeNode>(buildFactorGraph());
    return sumProduct.logNormalization();
  }
  
  /**
   * Next, we will use the matrix exponential and a discrete factor graph 
   * to compute the likelihood of the data under a certain tree.
   * 
   * First, derive the form of the factor graph that is needed to do this calculation.
   * See lecture 10 for a refresher on the model.
   * 
   * Look at the documentation in DiscreteFactorGraph for help on how to 
   * create the factor graph.
   * 
   * Look also at EdgeSorter to orient the edges of the tree (recall that you 
   * can root the tree arbitrarily).
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
  
  /**
   * Creates a tree likelihood factor from the sequences in the provided file.
   * 
   * The format of the file is: 
   *  - one line per taxon (species/leaf in the tree)
   *  - each line should begin with the name of the taxon, a space, and then the DNA sequence
   *  
   * See for example the file primates.data in this repository for an example.
   * 
   * The default rate matrix used is the Juke Cantor matrix, giving the same rate
   * to all transitions. 
   * 
   * @param inputFile The file containing the sequence data.
   */
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
  
  private TreeLikelihood(UnrootedTree tree, R rateMatrix,
      Map<TreeNode, UnaryFactor<TreeNode>> observationFactors, int nSites)
  {
    this.tree = tree;
    this.rateMatrix = rateMatrix;
    this.observationFactors = observationFactors;
    this.nSites = nSites;
  }
}
