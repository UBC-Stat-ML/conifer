package conifer.mmpp;

import bayonet.distributions.Exponential;
import bayonet.distributions.Uniform;
import bayonet.math.EJMLUtils;
import briefj.BriefIO;
import briefj.run.Results;
import conifer.ctmc.*;
import org.apache.commons.math3.util.Pair;
import org.ejml.simple.SimpleMatrix;
import org.junit.Test;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Created by crystal on 2016-04-20.
 */
public class TestForwardAndEndPointSamplersMMPP
{

    public static void main(String [] args)
    {
        new TestForwardAndEndPointSamplersMMPP().testUsingFixedRateMtxAndIntensity();
        new TestForwardAndEndPointSamplersMMPP().testUsingRandomRateMtxAndIntensity();
    }


    @Test
    public void testUsingRandomRateMtxAndIntensity(){

        //Generate a random GTR matrix and a random intensity vector, then test
        Random rand = new Random(1);
        int nStates = 5;
        SimpleRateMatrix gtr = RateMatrices.randomGTR(rand, nStates);
        System.out.println(gtr.toString());
        System.out.println();


        //Generate intensity, a vector of five elements
        double [] intensities = new double[nStates];
        double [] seeds = {1, 2, 3, 4, 5};
        for(int i=0; i<nStates; i++)
        {
            intensities[i] = Exponential.generate(rand, seeds[i]);
        }

        System.out.println(Arrays.toString(intensities));
        System.out.println();

        final double T = 3.0;
        PathStatistics pathStat1 = new PathStatistics((gtr.nStates()+1));
        PathStatistics pathStat2 = new PathStatistics((gtr.nStates()+1));

        MMPP mmpp = new EigenCTMCMMPP(gtr.getRateMatrix(), intensities);

        ForwardSamplerMMPP fwdSampler = new ForwardSamplerMMPP(mmpp);

        EndPointSamplerMMPP postSampler = new EndPointSamplerMMPP(mmpp);

        int nIters = 1000000;
        Pair<SimpleMatrix, SimpleMatrix> result = forwardAndBackwardPathStat(fwdSampler, postSampler, T, nIters,
                rand, pathStat1, pathStat2, gtr);

        SimpleMatrix m1 = result.getFirst();
        SimpleMatrix m2 = result.getSecond();

        System.out.println(EJMLUtils.toString(m1));
        System.out.println();

        System.out.println(EJMLUtils.toString(m2));
        EJMLUtils.checkIsClose(m1, m2, 0.05);

        System.out.println(pathStat1.toString());
        System.out.println();

        System.out.println(pathStat2.toString());
        System.out.println();

    }

    @Test
    public void testUsingFixedRateMtxAndIntensity()
    {

        SimpleRateMatrix k80 = RateMatrices.kimura1980();
        Random rand = new Random(1);
        final double T = 3.0;
        int nIters = 1000000;
        PathStatistics pathStat1 = new PathStatistics((k80.nStates()+1));
        PathStatistics pathStat2 = new PathStatistics((k80.nStates()+1));

        CTMC process = k80.getProcess();
        double[] intensities = {0.1, 0.2, 0.3, 0.4};
        MMPP mmpp = new EigenCTMCMMPP(k80.getRateMatrix(), intensities);

        ForwardSamplerMMPP fwdSampler = new ForwardSamplerMMPP(mmpp);

        EndPointSamplerMMPP postSampler = new EndPointSamplerMMPP(mmpp);


        Pair<SimpleMatrix, SimpleMatrix> result = forwardAndBackwardPathStat(fwdSampler, postSampler, T, nIters,
                rand, pathStat1, pathStat2, k80);

        SimpleMatrix m1 = result.getFirst();
        SimpleMatrix m2 = result.getSecond();
        System.out.println(EJMLUtils.toString(m1));
        System.out.println();

        System.out.println(EJMLUtils.toString(m2));
        EJMLUtils.checkIsClose(m1, m2, 0.05);

        System.out.println(pathStat1.toString());
        System.out.println();

        System.out.println(pathStat2.toString());
        System.out.println();


    }


    public Pair<SimpleMatrix, SimpleMatrix> forwardAndBackwardPathStat(ForwardSamplerMMPP
            fwdSampler, EndPointSamplerMMPP postSampler, double T, int nIters, Random rand,
            PathStatistics pathStat1, PathStatistics pathStat2,SimpleRateMatrix gtr){


        for (int i = 0; i < nIters; i++)
        {
            MMPPPath current = new MMPPPath();
            fwdSampler.sample(rand, T, pathStat1, current);

            // get the start and end state based on current path
            int startState = current.firstState();
            int endState = current.lastState();
            boolean lastEndPointPoissonEvent = false;
            if(endState == gtr.nStates())
                lastEndPointPoissonEvent = true;

            List<Integer> listOfStates = current.getPoissonEventState();
            List<Double> listOfTimes = current.getPoissonTimes();
            int numPaths= listOfStates.size();

            if(numPaths==0)
            {
                Path p2 = new Path();
                postSampler.sample(rand, startState, endState, T, pathStat2, p2, lastEndPointPoissonEvent);
            }else{
                // perform endpoint sampler on each segment of the paths
                int beginningState = startState;
                double startingTime =0;
                double time=0;

                for(int j=0; j<numPaths;j++)
                {
                    Path p2 = new Path();
                    int finishingStateOfPath = listOfStates.get(j);
                    double finishingTime = listOfTimes.get(j);
                    time = finishingTime - startingTime;
                    postSampler.sample(rand, beginningState, finishingStateOfPath, time, pathStat2, p2, true);
                    startingTime=finishingTime;
                    beginningState = finishingStateOfPath;
                }
                Path p2 = new Path();
                double lastTimes = listOfTimes.get((listOfTimes.size()-1));
                postSampler.sample(rand, beginningState, endState, T-lastTimes, pathStat2, p2, lastEndPointPoissonEvent);

            }
            System.out.println(i);

        }

        SimpleMatrix m1 = pathStat1.getCountsAsSimpleMatrix().scale(1.0/nIters);
        SimpleMatrix m2 = pathStat2.getCountsAsSimpleMatrix().scale(1.0/nIters);

        Pair<SimpleMatrix, SimpleMatrix> result = new Pair<SimpleMatrix, SimpleMatrix>(m1, m2);
        return result;
    }

}
