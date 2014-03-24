package conifer.factors;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.models.MultiCategorySubstitutionModel;
import blang.ForwardSampler;
import blang.MCMCAlgorithm;
import blang.MCMCRunner;
import blang.annotations.DefineFactor;



public class TestUnrootedTreeLikelihood extends MCMCRunner
{
  
  @DefineFactor
  public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel> likelihood = 
    UnrootedTreeLikelihood.createEmptyDefaultLikelihood(1, TopologyUtils.syntheticTaxaList(4));
  
  
  public static void main(String [] args)
  {
    TestUnrootedTreeLikelihood runner = new TestUnrootedTreeLikelihood();
    runner.factory.setCheckAllNodesCoveredByMCMCMoves(false);
    
    Map<String,Double> prs = Maps.newHashMap();
    for (int i = 0; i < 10000; i++)
    {
      MCMCAlgorithm mcmc = runner.buildMCMCAlgorithm();
      ForwardSampler forwardSampler = new ForwardSampler(mcmc.model);
      forwardSampler.simulate(mcmc.options.random);
      String current = runner.likelihood.observations.toString();
      double pr = Math.exp( runner.likelihood.logDensity() );
      prs.put(current, pr);
  //    System.out.println(runner.likelihood.observations);
    }
    System.out.println(prs);
    double sum = 0.0;
    for (double value : prs.values())
      sum += value;
    System.out.println(sum);
  }
}
