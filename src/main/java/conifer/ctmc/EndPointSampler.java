package conifer.ctmc;

import java.util.List;
import java.util.Random;

import org.ejml.simple.SimpleMatrix;

import bayonet.distributions.Exponential;
import bayonet.distributions.Multinomial;
import bayonet.math.SpecialFunctions;

import com.google.common.collect.Lists;


/**
 * WARNING: Not thread-safe.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 * The objective of this code is to sample the full path given the initial states and the final states.
 * The number of transitions and sojourn time is unknown and needs to be sampled.
 * The theory behind this code can be found at 
 * http://www.stat.ubc.ca/~bouchard/courses/stat547-sp2013-14/lecture/2014/02/05/lecture10.html
 *
 */
public class EndPointSampler
{
  private static final int MAX_N_TRANSITION = 1000000;
  private final CTMC ctmc;
  private final SimpleMatrix uniformizedTransition;
  public final double maxDepartureRate;
  private final List<SimpleMatrix> cache;
  private double [] sojournWorkArray = new double[10];
  private final double [] transitionWorkArray;
  
  public EndPointSampler(CTMC ctmc)
  {
    this.ctmc = ctmc;
    this.maxDepartureRate = maxDepartureRate(ctmc.getRateMatrix());
    this.uniformizedTransition = new SimpleMatrix(uniformizedTransition(ctmc.getRateMatrix(), maxDepartureRate));
    this.cache = initCache();
    this.transitionWorkArray = new double[ctmc.getRateMatrix().length];
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
  public void sample(Random rand, int startPoint, int endPoint, double T, PathStatistics statistics, Path path)
  {
    if (path != null && !path.isEmpty() && path.lastState() != startPoint)
      throw new RuntimeException("Incompatible extension of the provided path");
    int nTransitions = sampleNTransitions(rand, startPoint, endPoint, T);
    generatePath(rand, startPoint, endPoint, T, nTransitions, path, statistics);
  }
  
  public void sample(Random rand, int startPoint, int endPoint, double T, PathStatistics statistics)
  {
    sample(rand, startPoint, endPoint, T, statistics, null);
  }
  
  public int cacheSize()
  {
    return cache.size();
  }
/**
 * This unifomizedTransition is used to create matrix "B" in the Alex's lecture notes 
 * @param rateMatrix
 * @param mu  can be understood as the omega in the lecture notes
 * @return
 */
  private static double[][] uniformizedTransition(double [][] rateMatrix, double mu)
  {
    int nStates = rateMatrix.length;
    double [][] result = new double[nStates][nStates];
    for (int i = 0; i < nStates; i++)
      for (int j = 0; j < nStates; j++)
        result[i][j] = (i == j ? 1.0 : 0.0) + rateMatrix[i][j] / mu;
    return result;
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
    double[][] transitionMarginal = ctmc.marginalTransitionProbability(T);
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
        Math.log(getUniformizedTransitionPower(nTransition).get(startPoint, endPoint));
      final double logDenom = SpecialFunctions.logFactorial(nTransition);
      final double current = Math.exp(logNum - logDenom);
      sum += current;
      if (sum >= uniform)
        return nTransition;
    }
    throw new RuntimeException("Max number of transitions exceeded " + MAX_N_TRANSITION);
  }
  
  private void generatePath(Random rand, int startPoint, int endPoint, double T, int nTransitions, Path resultPath, PathStatistics stat)
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
          getUniformizedTransitionPower(nTransitions - transitionIndex - 1).get(candidateState, endPoint);
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
  
  private SimpleMatrix getUniformizedTransitionPower(int power)
  {
    ensureCache(power);
    return cache.get(power);
  }

  private List<SimpleMatrix> initCache()
  {
    List<SimpleMatrix> result = Lists.newArrayList();
    result.add(SimpleMatrix.identity(uniformizedTransition.numCols()));
    return result;
  }

  private void ensureCache(int power)
  {
    int maxPowerInCache = cache.size() - 1;
    for (int curPower = maxPowerInCache + 1; curPower <= power; curPower++)
      cache.add(uniformizedTransition.mult(cache.get(curPower-1)));
  }
}
