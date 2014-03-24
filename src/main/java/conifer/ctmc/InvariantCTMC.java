package conifer.ctmc;

import org.ejml.simple.SimpleMatrix;

import bayonet.math.EJMLUtils;



public class InvariantCTMC implements CTMC
{
  private final double [] stationary;
  
  public InvariantCTMC(double [] stationary)
  {
    this.stationary = stationary;
  }

  @Override
  public double[][] marginalTransitionProbability(double branchLength)
  {
    return identity();
  }
  
  private double[][] _identity = null;
  private double[][] identity()
  {
    if (_identity == null)
      _identity = EJMLUtils.copyMatrixToArray(SimpleMatrix.identity(nStates()));
    return _identity;
  }

  private int nStates()
  {
    return stationary.length;
  }

  @Override
  public double[] stationaryDistribution()
  {
    return stationary;
  }

  @Override
  public double[][] getRateMatrix()
  {
    return new double[stationary.length][stationary.length];
  }
}
