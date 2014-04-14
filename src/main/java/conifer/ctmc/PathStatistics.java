package conifer.ctmc;

import java.util.Arrays;

import org.ejml.data.DenseMatrix64F;
import org.ejml.simple.SimpleMatrix;

import bayonet.math.EJMLUtils;



public class PathStatistics
  {
    private final double[] initialCounts;
    private final DenseMatrix64F counts;
    
    public SimpleMatrix getCountsAsSimpleMatrix()
    {
      return new SimpleMatrix(counts);
    }
    
    public PathStatistics(int nStates)
    {
      this.initialCounts = new double[nStates];
      this.counts = new DenseMatrix64F(nStates, nStates);
    }
    
    public double getSojournTime(int state)
    {
      return counts.get(state, state);
    }
    
    public void addSojournTime(int currentPoint, double time)
    {
      counts.add(currentPoint, currentPoint, time);
    }
    
    public double getTransitionCount(int currentState, int nextState)
    {
      if (currentState == nextState)
        throw new RuntimeException();
      return counts.get(currentState, nextState);
    }

    public void addTransition(int currentState, int nextState)
    {
      if (currentState != nextState)
        counts.add(currentState, nextState, 1.0);
    }
    
    public double getInitialCount(int state)
    {
      return initialCounts[state];
    }
    
    public void addInitial(int state)
    {
      initialCounts[state]++;
    }
    
    @Override
    public String toString()
    {
      return "initialCounts: " + Arrays.toString(initialCounts) + "\n" +
      		"transitionsAndSojourns:\n" + EJMLUtils.toString( new SimpleMatrix(counts));
    }
    
  }