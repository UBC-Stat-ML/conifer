package conifer;

import org.junit.Test;

import bayonet.distributions.Exponential.RateParameterization;
import blang.MCMCAlgorithm;
import blang.MCMCRunner;
import blang.annotations.DefineFactor;
import blang.validation.CheckStationarity;
import conifer.ctmc.JukeCantorRateMatrix;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;


/**
 * Test the phylogenetic MCMC moves on a simple tree model.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class TestSimplePhyloModel extends MCMCRunner
{
  @DefineFactor
  UnrootedTreeLikelihood<JukeCantorRateMatrix> treeLikelihood = UnrootedTreeLikelihood.createEmptySyntheticNucleotideLikelihood(2, TopologyUtils.syntheticTaxaList(5));
  
  @DefineFactor
  NonClockTreePrior<RateParameterization> treePrior = NonClockTreePrior.on(treeLikelihood.tree);

  @Test
  public void checkStationarity()
  {
    this.factory.setCheckAllNodesCoveredByMCMCMoves(false);
    MCMCAlgorithm algo = buildMCMCAlgorithm();
    

    System.out.println(algo);
    
    algo.options.nMCMCSweeps = 10;
    
    // Actual code for setting up the test itself
    CheckStationarity check = new CheckStationarity();
    check.setShowSampleSummaryStats(true);
    System.out.println("Summary statistics of the samples");
    System.out.println("---------------------------------");
    
    // Here: 1000 is the number of test iterations (different than MCMC sweeps, see below)
    //       0.05 is a p-value threshold
    check.check(algo, 10000, 0.05);
    
//    next: 
//      - debug to check why TreeLenPro stat are NaN
//      - investigate one sampler at the time
//      - Later: notes simplify the exclude/move mechanism, Variable class
    
  }
}
