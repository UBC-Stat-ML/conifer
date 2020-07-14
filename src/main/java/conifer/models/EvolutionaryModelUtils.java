package conifer.models;

import java.util.List;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.marginal.FactorGraph;
import bayonet.marginal.UnaryFactor;
import bayonet.marginal.algo.SumProduct;

import com.google.common.collect.Lists;

import conifer.EvolutionaryModel;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.io.TreeObservations;



public class EvolutionaryModelUtils
{
  /**
   * Build a factor graph for the given tree according to a model and observation.
   * 
   * @param evolutionaryModel
   * @param tree
   * @param root
   * @param observations
   * @return
   */
  public static List<FactorGraph<TreeNode>> buildFactorGraphs( 
      EvolutionaryModel evolutionaryModel,
      UnrootedTree tree,
      TreeNode root,
      TreeObservations observations)
  {
    return buildFactorGraphs(evolutionaryModel, tree, root, observations, true);
  }
  public static List<FactorGraph<TreeNode>> buildFactorGraphs( 
      EvolutionaryModel evolutionaryModel,
      UnrootedTree tree,
      TreeNode root,
      TreeObservations observations,
      boolean useInitialDistribution)
  {
    List<FactorGraph<TreeNode>> result = Lists.newArrayList();
    List<Pair<TreeNode,TreeNode>> orientedEdges = tree.getRootedEdges(root); 
    
    for (int i = 0; i < evolutionaryModel.nFactorGraphs(); i++)
    {
      FactorGraph<TreeNode> factorGraph = evolutionaryModel.newFactorGraph(tree.getTopology());
      LikelihoodFactoryContext context = new LikelihoodFactoryContext(factorGraph, tree, observations, i, 1.0);
      
      // add observations
      if (observations != null)
        for (TreeNode observedNode : observations.getObservedTreeNodes())
          if (tree.getTopology().vertexSet().contains(observedNode)) // this check is useful when the data contains nodes not in the tree (see for example SRPMove for a use case)
            evolutionaryModel.buildObservation(observedNode, context);
        
      // initial distribution
      if (useInitialDistribution)
        evolutionaryModel.buildInitialDistribution(root, context);
      
      // add transitions
      for (Pair<TreeNode,TreeNode> edge : orientedEdges)
        evolutionaryModel.buildTransition(edge.getLeft(), edge.getRight(), context);
      
      result.add(factorGraph);
    }
    
    return result;
  }
  
  public static List<SumProduct<TreeNode>> getSumProductsFromFactorGraphs(List<FactorGraph<TreeNode>> factorGraphs, TreeNode arbitraryRoot)
  {
    List<SumProduct<TreeNode>> sumProds = Lists.newArrayList();
    for (FactorGraph<TreeNode> factorGraph : factorGraphs)
      sumProds.add(new SumProduct<TreeNode>(factorGraph));
    return sumProds;
  }
  
  public static List<UnaryFactor<TreeNode>> getRootMarginalsFromFactorGraphs(List<FactorGraph<TreeNode>> factorGraphs, TreeNode arbitraryRoot)
  {
    return getRootMarginals(getSumProductsFromFactorGraphs(factorGraphs, arbitraryRoot), arbitraryRoot);
  }
  public static List<UnaryFactor<TreeNode>> getRootMarginals(List<SumProduct<TreeNode>> sumProds, TreeNode arbitraryRoot)
  {
    List<UnaryFactor<TreeNode>> result = Lists.newArrayList();
    for (SumProduct<TreeNode> sumProd : sumProds)
      result.add(sumProd.computeMarginal(arbitraryRoot));
    return result;
  }
  
}
