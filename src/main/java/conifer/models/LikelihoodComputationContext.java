package conifer.models;

import java.util.List;

import bayonet.marginal.FactorGraph;
import bayonet.marginal.UnaryFactor;
import bayonet.marginal.algo.SumProduct;

import com.google.common.collect.Lists;

import conifer.TreeNode;

/**
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class LikelihoodComputationContext
{
  public LikelihoodComputationContext(List<FactorGraph<TreeNode>> factorGraphs, TreeNode arbitraryRoot)
  {
    this.factorGraphs = factorGraphs;
    this.arbitraryRoot = arbitraryRoot;
    this.rootMarginals = null;
  }
  
  public LikelihoodComputationContext(List<UnaryFactor<TreeNode>> rootMarginals)
  {
    this.rootMarginals = rootMarginals;
    this.arbitraryRoot = null;
    this.factorGraphs = null;
  }
  
  // standard way:
  private final List<FactorGraph<TreeNode>> factorGraphs;
  private final TreeNode arbitraryRoot;
  
  // or more directly(useful when efficiently visiting neighborhood systems)
  private final List<UnaryFactor<TreeNode>> rootMarginals;
  
  private boolean useDirectSpec()
  {
    boolean direct = rootMarginals != null;
    boolean check = factorGraphs != null && arbitraryRoot != null;
    if ((direct && check) || (!direct && ! check))
      throw new RuntimeException();
    return direct;
  }
  
  public List<UnaryFactor<TreeNode>> getRootMarginals()
  {
    if (useDirectSpec())
      return rootMarginals;
    else
      return EvolutionaryModelUtils.getRootMarginalsFromFactorGraphs(factorGraphs, arbitraryRoot);
  }
  
  
}