package conifer;

import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import org.junit.Test;

import bayonet.distributions.Exponential.RateParameterization;
import blang.ForwardSampler;
import blang.MCMCAlgorithm;
import blang.MCMCRunner;
import blang.annotations.DefineFactor;
import blang.validation.CheckStationarity;
import briefj.opt.Option;
import conifer.ctmc.cnv.CopyNumberMixture;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.io.CopyNumberTreeObservation;
import conifer.models.CNMultiCategorySubstitutionModel;
import conifer.moves.AllBranchesScaling;
import conifer.moves.SPRMove;


/**
 * Test the phylogenetic MCMC moves on a simple tree model.
 * 
 * @author Sean Jewell (jewellsean@gmail.com)
 *
 */
public class TestCNPhyloModel extends MCMCRunner
{
    @Option
    public final int nSites = 100;
    @Option
    public final int nTaxa = 3;

    Set<TreeNode> leaves = new HashSet<TreeNode>(TopologyUtils.makeLeaves(nTaxa, "t_"));

    @DefineFactor
    public final UnrootedTreeLikelihood<CNMultiCategorySubstitutionModel<CopyNumberMixture>> likelihood = UnrootedTreeLikelihood.createEmptyCN(nSites, leaves);

    @DefineFactor
    NonClockTreePrior<RateParameterization> treePrior = NonClockTreePrior.on(likelihood.tree);

    @Test
    public void checkStationarity()
    {
        this.factory.mcmcOptions.random = new Random(10001);
        this.factory.mcmcOptions.CODA = false;
        this.factory.setCheckAllNodesCoveredByMCMCMoves(false);
        this.factory.excludeNodeMove(AllBranchesScaling.class);
        this.factory.excludeNodeMove(SPRMove.class);

        MCMCAlgorithm algo = buildMCMCAlgorithm();
        
        // TODO: check that need data for part of the testing? 
        ForwardSampler f = new ForwardSampler(algo.model);
        f.simulate(algo.options.random);
        ((CopyNumberTreeObservation) likelihood.observations).init();
        
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
