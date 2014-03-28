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
    public void addSojournTime(int currentPoint, double time)
    {
      counts.add(currentPoint, currentPoint, time);
    }

    public void addTransition(int currentState, int nextState)
    {
      if (currentState != nextState)
        counts.add(currentState, nextState, 1.0);
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