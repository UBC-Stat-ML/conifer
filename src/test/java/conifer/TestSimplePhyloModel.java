package conifer;

import java.util.Random;

import org.junit.Test;

import bayonet.distributions.Exponential.RateParameterization;
import blang.MCMCAlgorithm;
import blang.MCMCRunner;
import blang.annotations.DefineFactor;
import blang.mcmc.MHMove;
import blang.validation.CheckStationarity;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.DiscreteGammaMixture;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.moves.AllBranchesScaling;
import conifer.moves.SPRMove;
import conifer.moves.SingleBranchScaling;


/**
 * Test the phylogenetic MCMC moves on a simple tree model.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class TestSimplePhyloModel extends MCMCRunner
{

  @DefineFactor
  public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<DiscreteGammaMixture>> likelihood = 
    UnrootedTreeLikelihood.createEmpty(1, TopologyUtils.syntheticTaxaList(4));
  
//  UnrootedTree tree = NonClockTreePrior.generateWithExponentialDistributedBranchLengths(new Random(1), 3);
  
  @DefineFactor
  NonClockTreePrior<RateParameterization> treePrior = NonClockTreePrior.on(likelihood.tree);

  @Test
  public void checkStationarity()
  {
    this.factory.mcmcOptions.random = new Random(10001);
    this.factory.setCheckAllNodesCoveredByMCMCMoves(false);
    
    this.factory.excludeNodeMove(AllBranchesScaling.class);
//    this.factory.excludeNodeMove(SPRMove.class);
    
    MCMCAlgorithm algo = buildMCMCAlgorithm();

    System.out.println(algo);
    
    algo.options.nMCMCSweeps = 20;
    
    // Actual code for setting up the test itself
    CheckStationarity check = new CheckStationarity();
    check.setShowSampleSummaryStats(true);
    System.out.println("Summary statistics of the samples");
    System.out.println("---------------------------------");
    
    // Here: 1000 is the number of test iterations (different than MCMC sweeps, see below)
    //       0.05 is a p-value threshold
    check.check(algo, 1000, 0.05);
    
  }
}
