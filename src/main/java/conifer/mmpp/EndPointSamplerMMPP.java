package conifer.mmpp;

import bayonet.distributions.Exponential;
import bayonet.distributions.Multinomial;
import bayonet.math.SpecialFunctions;
import com.google.common.collect.Lists;
import conifer.ctmc.CTMC;
import conifer.ctmc.Path;
import conifer.ctmc.PathStatistics;
import conifer.ctmc.RateMatrixUtils;
import org.jblas.DoubleMatrix;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Created by crystal on 2016-04-20.
 */
public class EndPointSamplerMMPP
{
    private static final int MAX_N_TRANSITION = 1000000;

    private final MMPP mmpp;

    private final DoubleMatrix uniformizedTransition;
    public final double maxDepartureRate;
    private final List<DoubleMatrix> cache;
    private double [] sojournWorkArray = new double[10];
    private final double [] transitionWorkArray;

    public static boolean cached;

    public EndPointSamplerMMPP(MMPP mmpp)
    {
        this.mmpp = mmpp;
        this.maxDepartureRate = maxDepartureRate(mmpp.getVirtualRateMatrix());
        this.uniformizedTransition = new DoubleMatrix(uniformizedTransition(mmpp.getRateMatrix(), maxDepartureRate, mmpp.getIntensities()));
        this.cache = initCache();
        this.transitionWorkArray = new double[mmpp.getRateMatrix().length];
    }

    /**
     * Note: adds transitions and sojourns into statistics, but not initial counts.
     *
     * @param rand
     * @param startPoint
     * @param endPoint
     * @param T
     * @param statistics
     * @param path
     */
    public void sample(Random rand, int startPoint, int endPoint, double T, PathStatistics statistics, Path path, boolean endPointPoissonEventOfMMPP)
    {
        if (path != null && !path.isEmpty() && path.lastState() != startPoint)
            throw new RuntimeException("Incompatible extension of the provided path");
        int nTransitions = sampleNTransitions(rand, startPoint, endPoint, T);
        generatePath(rand, startPoint, endPoint, T, nTransitions, path, statistics, endPointPoissonEventOfMMPP);
    }

    public void sample(Random rand, int startPoint, int endPoint, double T, PathStatistics statistics, boolean endPointPoissonEventOfMMPP)
    {
        sample(rand, startPoint, endPoint, T, statistics, null, endPointPoissonEventOfMMPP);
    }

    public int cacheSize()
    {
        return cache.size();
    }



    /**
     * This unifomizedTransition is used to create the transition probability matrix for uniformization method of MMPPs
     * @param rateMatrix
     * @param mu  can be understood as the intensity, which must be bigger or equal to the absolute value of the diagonal
     *            elements of the virtual rate matrix, Q-Lambda. Lambda is the diagonal matrix of the intensities for
     *            Poisson events of MMPPs
     * @return
     */
    private static double[][] uniformizedTransition(double [][] rateMatrix, double mu, double [] intensities)
    {
        int nStates = rateMatrix[0].length;
        double [][] QStar = RateMatrixUtils.rateMatrixMinusDiagIntensity(rateMatrix, intensities);
        double [][] result = new double[nStates][nStates];
        for (int i = 0; i < nStates; i++)
            for (int j = 0; j < nStates; j++)
                result[i][j] = (i == j ? 1.0 : 0.0) + QStar[i][j] / mu;
        return result;
    }



    private static double [][] uniformizedTransition(MMPP mmpp)
    {
        double [][] rateMatrix= mmpp.getRateMatrix();
        double [] intensities=mmpp.getIntensities();
        double [][] rateMatrixMinusIntensity = RateMatrixUtils.rateMatrixMinusDiagIntensity(rateMatrix, intensities);
        double mu= maxDepartureRate(rateMatrixMinusIntensity);
        return uniformizedTransition(rateMatrix, mu, intensities);

    }

    private static double [][] uniformizedTransition(double [][]rateMatrix, double []intensities)
    {
        double [][] rateMatrixMinusIntensity = RateMatrixUtils.rateMatrixMinusDiagIntensity(rateMatrix, intensities);
        double mu= maxDepartureRate(rateMatrixMinusIntensity);
        return uniformizedTransition(rateMatrix, mu, intensities);
    }

    public static void main(String [] args){
        double [][] rateMatrix={{-1, 0.2, 0.6, 0.2}, {0.2, -1, 0.2, 0.6}, {0.6, 0.2, -1, 0.2}, {0.2, 0.6, 0.2, -1}};
        double [] intensities={0.1, 0.2, 0.3, 0.4};
        MMPP mmpp = new EigenCTMCMMPP(rateMatrix, intensities);
        double [][] rateMatrixMinusDiagIntensity = RateMatrixUtils.rateMatrixMinusDiagIntensity(rateMatrix, intensities);
        System.out.println(Arrays.deepToString(rateMatrix));
        System.out.println(Arrays.deepToString(rateMatrixMinusDiagIntensity));
        double [][] result = uniformizedTransition(rateMatrix, intensities);
        System.out.println(Arrays.deepToString(result));
        System.out.println(Arrays.deepToString(mmpp.marginalTransitionProbability(rateMatrix, intensities, 1.0)));
    }



    /**
     * This code will return the biggest diagonal element in the rate matrix
     * @param rateMatrix
     * @return
     */

    private static double maxDepartureRate(double [][] rateMatrix)
    {
        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < rateMatrix.length; i++)
        {
            final double current = Math.abs(rateMatrix[i][i]);
            if (current > max)
                max = current;
        }
        return max;
    }

    /**
     * The objective of this code is to generate the number of transitions between different states
     *  http://www.stat.ubc.ca/~bouchard/courses/stat547-sp2013-14/lecture/2014/02/05/lecture10.html is a good reference
     *  To understand this code, check Equation 17-20 in lecture 10
     *  To get the marginal distribution, we marginalize the product of the poission part and the factor graph part
     *  In the code below, we separate the terms related to "n"- number of transition steps and those not so that we have logConstant and logMuT, not
     *  related to "n"
     * @param rand
     * @param startPoint
     * @param endPoint
     * @param T
     * @return
     */

    public int sampleNTransitions(Random rand, int startPoint, int endPoint, double T)
    {
        double[][] transitionMarginal = mmpp.marginalTransitionProbability(mmpp.getRateMatrix(), mmpp.getIntensities(), T);
        if (transitionMarginal[startPoint][endPoint] == 0)
            throw new RuntimeException();
        if (transitionMarginal[startPoint][endPoint] == 1.0)
            return 0;
        final double uniform = rand.nextDouble();
        double sum = 0.0;
        final double logConstant = (-maxDepartureRate * T) - Math.log( transitionMarginal[startPoint][endPoint] );
        final double logMuT = Math.log(maxDepartureRate * T);
        for (int nTransition = 0; nTransition < MAX_N_TRANSITION; nTransition++)
        {
            final double logNum =
                    logConstant +
                            nTransition * logMuT +
                            Math.log(getUniformizedTransitionPower(nTransition,cached).get(startPoint, endPoint));
            final double logDenom = SpecialFunctions.logFactorial(nTransition);
            final double current = Math.exp(logNum - logDenom);
            sum += current;
            if (sum >= uniform)
                return nTransition;
        }
        throw new RuntimeException("Max number of transitions exceeded " + MAX_N_TRANSITION);
    }

    private void generatePath(Random rand, int startPoint, int endPoint, double T, int nTransitions, Path resultPath, PathStatistics stat, boolean endPointPoissonEventOfMMPP)
    {
        double [] sojournTimes = generateSojournTimes(rand, nTransitions, T);
        int currentPoint = startPoint;
        if (resultPath != null)
            resultPath.addSegment(currentPoint, sojournTimes[0]);
        if (stat != null)
            stat.addSojournTime(currentPoint, sojournTimes[0]);
        for (int transitionIndex = 0; transitionIndex < nTransitions; transitionIndex++)
        {
            // compute transition probabilities
            for (int candidateState = 0; candidateState < transitionWorkArray.length; candidateState++)
                transitionWorkArray[candidateState] =
                        uniformizedTransition.get(currentPoint, candidateState) *
                                getUniformizedTransitionPower(nTransitions - transitionIndex - 1,cached).get(candidateState, endPoint);
            Multinomial.normalize(transitionWorkArray);
            int nextState = Multinomial.sampleMultinomial(rand, transitionWorkArray);
            if (resultPath != null)
                resultPath.addSegment(nextState, sojournTimes[transitionIndex+1]);
            if (stat != null)
            {
                stat.addSojournTime(nextState, sojournTimes[transitionIndex+1]);
                if (currentPoint != nextState)
                    stat.addTransition(currentPoint, nextState);
            }
            currentPoint = nextState;
        }
        if(endPointPoissonEventOfMMPP){
            stat.addTransition(endPoint, mmpp.getRateMatrix()[0].length);
        }
    }


    /**
     * Warning: for efficiency reason, the returned array may be longer than
     * the actual number of relevant cells
     * @param rand
     * @param nTransitions
     * @param T
     * @return
     */
    private double[] generateSojournTimes(Random rand, int nTransitions, double T)
    {
        final int nTimes = nTransitions + 1;
        double [] result = getWorkArray(nTimes);
        double sum = 0.0;
        for (int i = 0; i < nTimes; i++)
        {
            final double cur = Exponential.generate(rand, 1.0);
            sum += cur;
            result[i] = cur;
        }
        for (int i = 0; i < nTimes; i++)
            result[i] = T * result[i] / sum;
        return result;
    }

    private double[] getWorkArray(int minLen)
    {
        if (sojournWorkArray.length < minLen)
            sojournWorkArray = new double[minLen * 2];
        return sojournWorkArray;
    }

    private DoubleMatrix getUniformizedTransitionPower(int power)
    {
        ensureCache(power);
        return cache.get(power);
    }

    private DoubleMatrix getUniformizedTransitionPower(int power, boolean cached)
    {
        if(cached){
            return getUniformizedTransitionPower(power);
        }else{
            DoubleMatrix result = DoubleMatrix.eye(uniformizedTransition.columns);
            for(int i=0; i< power; i++)
            {
                result = result.mmul(uniformizedTransition);

            }
            return result;

        }

    }

    private List<DoubleMatrix> initCache()
    {
        List<DoubleMatrix> result = Lists.newArrayList();
        result.add(DoubleMatrix.eye(uniformizedTransition.columns));
        return result;
    }

    private void ensureCache(int power)
    {
        int maxPowerInCache = cache.size() - 1;
        for (int curPower = maxPowerInCache + 1; curPower <= power; curPower++)
        {
            DoubleMatrix tmpMatrix = cache.get(curPower-1);
            cache.add(tmpMatrix.mmul(uniformizedTransition));
        }

    }


}
