package conifer.moves;

import java.util.Map;
import java.util.Random;

import blang.mcmc.ConnectedFactor;
import blang.mcmc.MHProposalDistribution;
import blang.mcmc.SampledVariable;
import briefj.opt.Option;
import conifer.Parsimony;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.io.CopyNumberTreeObservation;
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
    
    @ConnectedFactor 
    
    @Override
    public Proposal propose(Random rand)
    {
        CopyNumberTreeObservation data = (CopyNumberTreeObservation) likelihood.observations;
               
        int nSites = data.nSites();
        for(int i = 0; i < nSites; i++)
        {
            Map<String, CNPair> emissions = data.getEmissionAtSite(i);
            int Mi = oneSiteGibbs(rand, emissions, likelihood);
            parsimony.getM().setIndex(i, Mi);
        }
        return new ProposalRealization();
    }
    
    
    private int oneSiteGibbs(Random rand, Map<String, CNPair> emissions, UnrootedTreeLikelihood<CNMultiCategorySubstitutionModel<RateMatrixMixture>> likelihood)
    {
        double[] loglikelihood; 
        int noLeaves = emissions.keySet().size(); 
        
        for (int i = 0; i < noLeaves; i++)
        {
            
        }
         
        
        loglikelihood[v] = likelihood.logDensity();
        
        return 0; 
    }
    
    
    private double[] conditionalEmissionLikelihood(int N, int rA, double xi, double s)
    {
        
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
