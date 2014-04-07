package conifer.ctmc.expfam;

import org.junit.Test;

import bayonet.distributions.Normal;
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
    UnrootedTreeLikelihood.createEmpty(1, TopologyUtils.syntheticTaxaList(2)).withExpFamMixture(ExpFamMixture.dnaGTR(2));
  
  @DefineFactor
  public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior =
    IIDRealVectorGenerativeFactor.iidNormalOn(likelihood.evolutionaryModel.rateMatrixMixture.parameters);
  
  
  @DefineFactor
  NonClockTreePrior<RateParameterization> treePrior = NonClockTreePrior.on(likelihood.tree);

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
    
  }
}
