package conifer.moves;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import blang.mcmc.ConnectedFactor;
import blang.mcmc.MHProposalDistribution;
import blang.mcmc.SampledVariable;
import briefj.Indexer;
import briefj.opt.Option;
import conifer.Parsimony;
import conifer.TreeNode;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.io.CopyNumberTreeObservation;
import conifer.io.Indexers;
import conifer.models.CNMultiCategorySubstitutionModel;
import conifer.models.CNPair;
import conifer.models.RateMatrixMixture;

/**
 * 
 * @author jewellsean
 *
 */
public class ParsimonySampler implements MHProposalDistribution
{
    @Option public final double DELTA = 0.1;

    @SampledVariable Parsimony parsimony; 

    @ConnectedFactor UnrootedTreeLikelihood<CNMultiCategorySubstitutionModel<RateMatrixMixture>> likelihood;

//    @ConnectedFactor need for s variable 

    @Override
    public Proposal propose(Random rand)
    {
        CopyNumberTreeObservation data = (CopyNumberTreeObservation) likelihood.observations;
         
        
        int nSites = data.nSites();
        for(int i = 0; i < nSites; i++)
        {
            Map<String, CNPair> emissions = data.getEmissionAtSite(i);
            int Mi = oneSiteGibbs(rand, emissions, i);
            parsimony.getM().setIndex(i, Mi);
        }
        return new ProposalRealization();
    }


    private int oneSiteGibbs(Random rand, Map<String, CNPair> emissions, int site)
    {
        double[] loglikelihood; 
        int noLeaves = emissions.keySet().size(); 

        for (int v = 0; v < noLeaves; v++)
        {
            updateAllLeavesFixedSite(emissions, v, betaBinomialPrecision, site);
            loglikelihood[v] = likelihood.logDensity();
        }

        int mi = bayonet.distributions.Multinomial.sampleMultinomial(rand, loglikelihood);
        updateAllLeavesFixedSite(emissions, mi, betaBinomialPrecision, site);
        
        return mi; 
    }

    // for a fixed parsimony state, that is a fixed value of M update all of the leaves with new Y vectors
    private void updateAllLeavesFixedSite(Map<String, CNPair> emissions, int minLeaf, double betaBinomialPrecision, int site)
    {
        CopyNumberTreeObservation data = (CopyNumberTreeObservation) likelihood.observations;
        Map<String, Integer> leafMap = data.getLeafOrder();
        for (String leaf : emissions.keySet())
        {
            double[] logLik = conditionalEmissionLogLikelihood(emissions.get(leaf), betaBinomialPrecision, leafMap.get(leaf), minLeaf);
            data.setSite(TreeNode.withLabel(leaf), site, logLik);
        }
    }


    // for a single leaf update the vector of emission probabilities 
    // imposes the conditions provided through M on Y vector
    private double[] conditionalEmissionLogLikelihood(CNPair emissions, double betaBinomialPrecision, int leaf, int minLeaf)
    {
        Indexer<String> indexer = Indexers.copyNumberCTMCIndexer();
        int noStates = indexer.objectsList().size();
        double[] logLik = new double[noStates];

        for (int i = 0; i < noStates; i++)
        {
            int[] state = parseState(indexer.i2o(i));
            int b = state[2];
            if(leaf > minLeaf || (b == 0 && leaf < minLeaf) || (b == 1) && leaf == minLeaf)
            {
                logLik[i] = conditionalEmissionLogLikelihood(emissions.getN(), emissions.getrA(), constructXi(state[0], state[1]), betaBinomialPrecision);    
            }
        }
        return logLik;
    }

    private int[] parseState(String state)
    {
        int stateSize = 3; 
        int[] stateVec = new int[3];

        String[] p = state.split(",");
        for (int i = 0; i < stateSize; i++)
        {
            stateVec[i] = Integer.valueOf(p[i]).intValue();    
        }

        return stateVec;

    }

    private double conditionalEmissionLogLikelihood(int trials, int rA, double xi, double s)
    {
        double alpha = s * xi; 
        double beta = s * (1 - xi);
        return bayonet.distributions.BetaBinomial.logDensity(rA, alpha, beta, trials);
    }

    private double constructXi(int A, int a)
    {
        return (A + DELTA) / (A + a + DELTA);
    }

    private class ProposalRealization implements Proposal
    {
        @Override
        public double logProposalRatio() // not symmetric, but Gibbs? 
        {
            return 0;
        }

        // This is a Gibbs move
        @Override
        public void acceptReject(boolean accept){}

    }

}