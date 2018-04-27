package conifer.ctmc;

import java.util.Random;

import bayonet.distributions.Multinomial;
import conifer.Utils;




public class ForwardSampler
{
  private final CTMC ctmc;
  
  public ForwardSampler(CTMC ctmc)
  {
    this.ctmc = ctmc;
  }

  public void sample(Random rand, double T, PathStatistics statistics, Path path)
  {
    int start = Multinomial.sampleMultinomial(rand, ctmc.stationaryDistribution());
    sample(rand, T, start, statistics, path);
  }
  
  public void sample(Random rand, double T, int startState, PathStatistics statistics, Path path)
  {
    double [][] rateMatrix = ctmc.getRateMatrix();
    double [][] embeddedMarkovChain = RateMatrixUtils.getJumpProcess(rateMatrix);
    
    double totalTime = 0.0;
    int state = startState;
    while (totalTime < T)
    {
      // sample waiting time
      double currentRate = -rateMatrix[state][state];
      double time = currentRate == 0 ? // absorbing state
          T - totalTime :
          Utils.sampleExponential(rand, currentRate);
      
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
      int nextState = Multinomial.sampleMultinomial(rand, embeddedMarkovChain[state]);
      if (!finished && statistics != null)
        statistics.addTransition(state, nextState);
      
      state = nextState;
    }
  }
}
