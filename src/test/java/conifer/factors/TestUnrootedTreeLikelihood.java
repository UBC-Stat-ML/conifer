package conifer.factors;

import java.util.Random;

import org.junit.Test;

import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.validation.CheckDiscreteNormalization;
import conifer.TopologyUtils;
import conifer.models.MultiCategorySubstitutionModel;



public class TestUnrootedTreeLikelihood
{
  @DefineFactor
  public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel> likelihood = 
    UnrootedTreeLikelihood.createEmptyDefaultLikelihood(1, TopologyUtils.syntheticTaxaList(4));
  
  public static void main(String [] args)
  {
    TestUnrootedTreeLikelihood runner = new TestUnrootedTreeLikelihood();
    runner.runTest();
  }
  
  @Test
  public void runTest()
  {
    
    ProbabilityModel model = new ProbabilityModel(this);
    Random rand = new Random(1);
    CheckDiscreteNormalization.check(model, rand, 10000);
//    runner.factory.setCheckAllNodesCoveredByMCMCMoves(false);
//    MCMCAlgorithm mcmc = runner.buildMCMCAlgorithm();
//    ForwardSampler forwardSampler = new ForwardSampler(mcmc.model);
//    Map<String,Double> prs = Maps.newHashMap();
//    for (int i = 0; i < 10000; i++)
//    {
//      forwardSampler.simulate(mcmc.options.random);
//      String current = runner.likelihood.observations.toString();
//      double pr = Math.exp( runner.likelihood.logDensity() );
//      prs.put(current, pr);
//  //    System.out.println(runner.likelihood.observations);
//    }
//    System.out.println(prs);
//    double sum = 0.0;
//    for (double value : prs.values())
//      sum += value;
//    System.out.println(sum);
  }
}
