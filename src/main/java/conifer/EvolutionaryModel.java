package conifer;

import java.util.Random;

import org.jgrapht.UndirectedGraph;

import bayonet.marginal.FactorGraph;
import briefj.collections.UnorderedPair;
import conifer.io.TreeObservations;
import conifer.models.LikelihoodComputationContext;
import conifer.models.LikelihoodFactoryContext;



public interface EvolutionaryModel
{
  public int nFactorGraphs();
  public FactorGraph<TreeNode> newFactorGraph(UndirectedGraph<TreeNode, UnorderedPair<TreeNode, TreeNode>> undirectedGraph);
  
  public void buildObservation(TreeNode leaf, LikelihoodFactoryContext context);
  public void buildTransition(TreeNode parent, TreeNode children, LikelihoodFactoryContext context);
  public void buildInitialDistribution(TreeNode node, LikelihoodFactoryContext context);
  
  public double computeLogLikelihood(LikelihoodComputationContext context);

  /**
   * Generates observations at the leaves of the TreeObservation object.
   * 
   * @param rand
   * @param destination
   * @param tree
   * @param root
   */
  public void generateObservationsInPlace(Random rand, TreeObservations destination, UnrootedTree tree, TreeNode root);
}