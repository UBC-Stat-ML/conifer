package conifer.ctmc;

import org.ejml.data.DenseMatrix64F;
import org.ejml.simple.SimpleMatrix;

import bayonet.math.EJMLUtils;



public class PathStatistics
  {
    private final DenseMatrix64F counts;
    
    public SimpleMatrix getCountsAsSimpleMatrix()
    {
      return new SimpleMatrix(counts);
    }
    
    public PathStatistics(int nStates)
    {
      this.counts = new DenseMatrix64F(nStates, nStates);
    }
//    private Counter<Pair<Integer, Integer>> transitionCounts = new Counter<Pair<Integer,Integer>>();
//    private Counter<Integer> sojournTimes = new Counter<Integer>();
    public void addSojournTime(int currentPoint, double time)
    {
      counts.add(currentPoint, currentPoint, time);
//      sojournTimes.incrementCount(currentPoint, d);
    }

    public void addTransition(int currentState, int nextState)
    {
      if (currentState != nextState)
        counts.add(currentState, nextState, 1.0);
//        transitionCounts.incrementCount(Pair.of(currentState, nextState), 1.0);
    }
    
    @Override
    public String toString()
    {
      return EJMLUtils.toString(new SimpleMatrix(counts));
    }
    
  }