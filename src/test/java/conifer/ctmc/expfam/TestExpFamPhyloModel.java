package conifer.ctmc.expfam;

import org.junit.Before;
import org.junit.Test;

import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.MCMCAlgorithm;
import blang.MCMCRunner;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.validation.CheckStationarity;
import conifer.TopologyUtils;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;


/**
 * Test the phylogenetic MCMC moves on a simple tree model.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class TestExpFamPhyloModel extends MCMCRunner
{
  @DefineFactor
  public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = 
    UnrootedTreeLikelihood.createEmpty(1, TopologyUtils.syntheticTaxaList(2)).withExpFamMixture(ExpFamMixture.dnaGTR());
  
  @DefineFactor
  public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior =
    IIDRealVectorGenerativeFactor.iidNormalOn(likelihood.evolutionaryModel.rateMatrixMixture.parameters);
  
  @DefineFactor
  NonClockTreePrior<RateParameterization> treePrior = NonClockTreePrior.on(likelihood.tree);

  public static void main(String [] args)
  {
//    shouldUsePlots = true;
    new TestExpFamPhyloModel().checkStationarity();
  }
  
  public static boolean shouldUsePlots = false;

  
  @Test
  public void checkStationarity()
  {
    this.factory.setCheckAllNodesCoveredByMCMCMoves(false);
    MCMCAlgorithm algo = buildMCMCAlgorithm();
    algo.options.CODA = shouldUsePlots;

    System.out.println(algo);
    
    algo.options.nMCMCSweeps = 10;
    
    // Actual code for setting up the test itself
    CheckStationarity check = new CheckStationarity();
    check.setShowSampleSummaryStats(true);
    System.out.println("Summary statistics of the samples");
    System.out.println("---------------------------------");
    
    // Here: 1000 is the number of test iterations (different than MCMC sweeps, see below)
    //       0.05 is a p-value threshold
    check.check(algo, 1000, 0.05);
    
    /*
     * Note: for MCITERS=10000,NCAT=2 (*) the test fails (small diff, but p-values in the 0.01)
     * but for   MCITERS=100000, NCAT=1 the test works
     * 
     * Seems to be due to the maximization-initialization of HMC,
     * as (*) does not fail when modifying HMC init to not keep the 
     * max used when finding hyper-param values.
     */
  }
}
