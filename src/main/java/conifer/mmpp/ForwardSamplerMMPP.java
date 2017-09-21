package conifer.mmpp;

import bayonet.distributions.Exponential;
import bayonet.distributions.Multinomial;
import conifer.ctmc.CTMC;
import conifer.ctmc.Path;
import conifer.ctmc.PathStatistics;
import conifer.ctmc.RateMatrixUtils;

import java.util.Random;

/**
 * Created by crystal on 2016-04-19.
 */
public class ForwardSamplerMMPP {
    private final MMPP mmpp;

    public ForwardSamplerMMPP(MMPP mmpp)
    {
        this.mmpp = mmpp;
    }

    public void sample(Random rand, double T, PathStatistics statistics, MMPPPath path)
    {
        int start = Multinomial.sampleMultinomial(rand, mmpp.stationaryDistribution());
        sample(rand, T, start, statistics, path);
    }

    public void sample(Random rand, double T, int startState, PathStatistics statistics, MMPPPath path)
    {
        double [][] rateMatrix = mmpp.getRateMatrix();
        double [] intensities = mmpp.getIntensities();
        // get the virtual rate matrix for MMPPs including the virtual state to represent Poisson Event
        double [][] virtualRateMatrix = RateMatrixUtils.fillVirtualRateMatrixForMMPPs(rateMatrix, intensities);
        double [][] embeddedMMPPs = RateMatrixUtils.getJumpProcessMMPP(virtualRateMatrix);

        double totalTime = 0.0;
        int state = startState;
        //startState can not be a Poisson event, so that it can not be the virtual state.
        int nCTMCState = intensities.length;
        //Since the index starts from zero, the index for the virtual state is nCTMCState;
        if(startState>=nCTMCState)
            throw new RuntimeException("The startState can not be a virtual state");


        while (totalTime < T)
        {
            // sample waiting time, use the diagonal elements from the virtual rate matrix of MMPPs
            double currentRate = -virtualRateMatrix[state][state];
            double time = currentRate == 0 ? // absorbing state
                    T - totalTime :
                    Exponential.generate(rand, currentRate);

            totalTime += time;
            boolean finished = false;
            if (totalTime > T)
            {
                time = time - (totalTime -T);
                finished = true;
            }

            if (statistics != null)
                statistics.addSojournTime(state, time);
            if (path != null)
                path.addSegment(state, time);

            // sample transition
            int nextState = Multinomial.sampleMultinomial(rand, embeddedMMPPs[state]);

            if (!finished && statistics != null)
                statistics.addTransition(state, nextState);

            // if the nextState is a virtual state, that is a Poisson event
            // Since we are using an embedded chain, if the next state is the same as the previous one,
            // this indicates it is a Poisson event.
            if(nextState==nCTMCState){
                nextState = state;
                if(totalTime<=T)
                {
                    path.addPoissonEventTimesOfMMPP(totalTime);
                    path.addPoissonEventStatesOfMMPP(state);
                }

            }
            state = nextState;
        }
    }
}
