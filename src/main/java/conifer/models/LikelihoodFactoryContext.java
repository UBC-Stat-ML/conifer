package conifer.models;

import java.util.Map;

import bayonet.marginal.DiscreteFactorGraph;
import bayonet.marginal.FactorGraph;

import com.google.common.collect.Maps;

import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.io.TreeObservations;

/**
 * Holds contextual information used by implementations of EvolutionaryModels
 * to construct one factor graph.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class LikelihoodFactoryContext
{
  private final FactorGraph<TreeNode> factorGraph;

  private final TreeObservations observations;
  private final int factorGraphIndex;
  
  private final double observationAnnealing;

  private final UnrootedTree tree;
  
  public LikelihoodFactoryContext(
      FactorGraph<TreeNode> factorGraph,
      UnrootedTree tree, TreeObservations observations, 
      int factorGraphIndex,
      double observationAnnealing)
  {
    this.factorGraph = factorGraph;
    this.observations = observations;
    this.factorGraphIndex = factorGraphIndex;
    this.tree = tree;
    this.observationAnnealing = observationAnnealing;
  }

  public TreeObservations getObservations()
  {
    return observations;
  }

  public DiscreteFactorGraph<TreeNode> getDiscreteFactorGraph()
  {
    return (DiscreteFactorGraph<TreeNode>) factorGraph;
  }
  
  @SuppressWarnings("rawtypes")
  private final Map cache = Maps.newHashMap();

  @SuppressWarnings("rawtypes")
  public Map getCache()
  {
    return cache;
  }

  public int getFactorGraphIndex()
  {
    return factorGraphIndex;
  }

  public double getBranchLength(TreeNode parent, TreeNode children)
  {
    return tree.getBranchLength(parent, children);
  }

  public double getObservationAnnealing() {
    return observationAnnealing;
  }

}