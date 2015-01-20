package conifer.moves;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import blang.mcmc.ConnectedFactor;
import blang.mcmc.MHProposalDistribution;
import blang.mcmc.SampledVariable;
import briefj.Indexer;
import briefj.opt.Option;
import conifer.TreeNode;
import conifer.ctmc.cnv.CopyNumberMixture;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.io.CopyNumberTreeObservation;
import conifer.io.Indexers;
import conifer.models.CNPair;
import conifer.models.MultiCategorySubstitutionModel;

/**
 * 
 * @author Sean Jewell (jewellsean@gmail.com)
 *
 */
public class CopyNumberTreeSampler implements MHProposalDistribution
{
    // TODO: this needs to be moved out
    @Option public static final double DELTA = 0.1;

    @SampledVariable CopyNumberTreeObservation tree; 

    @ConnectedFactor UnrootedTreeLikelihood<MultiCategorySubstitutionModel<CopyNumberMixture>> likelihood;

    private List<Integer> visitSiteInRandomOrder(Random rand, int nSites)
    {
        List<Integer> visitOrder = new ArrayList<Integer>();
        for (Integer i = 0; i < nSites; i++)
            visitOrder.add(i);
        Collections.shuffle(visitOrder, rand);
        return visitOrder;
    }

    @Override
    public Proposal propose(Random rand)
    {
        CopyNumberTreeObservation data = (CopyNumberTreeObservation) likelihood.observations;
        List<Integer> visitOrder = visitSiteInRandomOrder(rand, data.nSites());

        for(Integer i : visitOrder)
        {
            Map<String, Set<CNPair>> emissions = data.getEmissionAtSite(i.intValue());
            int Mi = oneSiteGibbs(rand, emissions, i.intValue());
            tree.parsimony.getM().setIndex(i.intValue(), Mi);
        }
        return new ProposalRealization();
    }


    private int oneSiteGibbs(Random rand, Map<String, Set<CNPair>> emissions, int site)
    {
        int noLeaves = emissions.keySet().size(); 
        double[] loglikelihood = new double[noLeaves];

        for (int v = 0; v < noLeaves; v++)
        {
            updateAllLeavesFixedSite(emissions, v, site);
            loglikelihood[v] = likelihood.logDensity();
        }
        double[] prob = bayonet.opt.DoubleArrays.pointwiseExp(loglikelihood);
        bayonet.distributions.Multinomial.normalize(prob);
        int mi = bayonet.distributions.Multinomial.sampleMultinomial(rand, prob);
        updateAllLeavesFixedSite(emissions, mi, site);

        return mi; 
    }

    // for a fixed parsimony state, that is a fixed value of M update all of the leaves with new Y vectors
    private void updateAllLeavesFixedSite(Map<String, Set<CNPair>> emissions, int minLeaf, int site)
    {
        CopyNumberTreeObservation data = (CopyNumberTreeObservation) likelihood.observations;
        Map<String, Integer> leafMap = data.getLeafOrder();
        for (String leaf : emissions.keySet())
        {
            double[] logLik = conditionalEmissionLogLikelihood(emissions.get(leaf), tree.betaBinomialprecision.getValue(), leafMap.get(leaf), minLeaf);
            data.setSite(TreeNode.withLabel(leaf), site, bayonet.opt.DoubleArrays.pointwiseExp(logLik));
        }
    }


    // for a single leaf update the vector of emission probabilities 
    // imposes the conditions provided through M on Y vector
    private double[] conditionalEmissionLogLikelihood(Set<CNPair> emissions, double betaBinomialPrecision, int leaf, int minLeaf)
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
                // loop over all elements in cluster 
                Iterator<CNPair> clusterElements = emissions.iterator();
                while(clusterElements.hasNext())
                {
                    CNPair currentElement = clusterElements.next();
                    logLik[i] += conditionalEmissionLogLikelihood(currentElement.getN(), currentElement.getrA(), constructXi(state[0], state[1]), betaBinomialPrecision);
                }

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

    // TODO: needs to be moved out of this class
    // TODO: validate reasonableness of this formualtion
    public static double constructXi(int A, int a)
    {
        double xi; 
        if (A == 0 && a == 0 )
            xi = 0.5;
        else if (A == 0)
            xi = DELTA / (a + DELTA);
        else if (a == 0)
            xi = A / (A + DELTA);
        else
            xi = ((double) A) / (A + a);
        return xi;
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
