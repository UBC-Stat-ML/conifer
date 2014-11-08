package conifer.moves;

import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import conifer.Parsimony;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.io.CopyNumberTreeObservation;
import conifer.models.CNMultiCategorySubstitutionModel;
import conifer.models.CNPair;
import conifer.models.CNSpecies;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.models.ParsimonyModel;
import conifer.models.RateMatrixMixture;
import blang.factors.Factor;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.MHProposalDistribution;
import blang.mcmc.SampledVariable;

/**
 * 
 * @author jewellsean
 *
 */
public class ParsimonySampler implements MHProposalDistribution
{
    @SampledVariable Parsimony parsimony; 
    
    @ConnectedFactor UnrootedTreeLikelihood<CNMultiCategorySubstitutionModel<RateMatrixMixture>> likelihood;
    
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
        
        
        
        
        loglikelihood[v] = likelihood.logDensity();
        
        return 0; 
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
