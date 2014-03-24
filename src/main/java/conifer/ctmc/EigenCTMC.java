package conifer.ctmc;

import org.ejml.data.Eigenpair;
import org.ejml.ops.EigenOps;
import org.ejml.simple.SimpleMatrix;

import tutorialj.Tutorial;
import bayonet.distributions.Multinomial;
import bayonet.math.EJMLUtils;
import bayonet.math.NumericalUtils;
import bayonet.math.EJMLUtils.SimpleEigenDecomposition;


/**
 * A continuous time Markov chain. The main functionalities consists
 * in computing a marginal transition probability and a stationary distibution
 * (see below).
 * 
 * This implementation is based on caching the eigendecomposition
 * of the provided rate matrix.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class EigenCTMC implements CTMC
{
  private final SimpleEigenDecomposition eigenDecomp;
  private final double [] stationaryDistribution;
  private final double [][] rates;

  /**
   * Note: if the RateMatrix is changed in place,
   * these changes will not be mirrored by this class.
   * 
   * It should be recreated each time a likelihood 
   * calculation is performed.
   * @param rateMatrix
   */
  public EigenCTMC(double [][] rates)
  {
    RateMatrixUtils.checkValidRateMatrix(rates);
    this.eigenDecomp = EJMLUtils.simpleEigenDecomposition(new SimpleMatrix(rates));
    this.stationaryDistribution = computeStationary();
    this.rates = rates;
  }
  
  /**
   * Compute and cache the stationary (called by the constructor)
   * @return
   */
  private double[] computeStationary()
  {
    SimpleMatrix marginal = new SimpleMatrix(marginalTransitionProbability(1.0));
    SimpleMatrix exponentiatedTranspose = marginal.transpose();
    Eigenpair eigenPair = EigenOps.computeEigenVector(exponentiatedTranspose.getMatrix(), 1.0);
    double [] result = eigenPair.vector.data;
    double sum = NumericalUtils.getNormalization(result);
    if (sum < 0.0)
      for (int i = 0; i < result.length; i++)
        result[i] *= -1;
    Multinomial.normalize(result);
    return result;
  }
  

  /**
   * 
   * Fill in ``marginalTransitionProbability()``.
   * 
   * Use the diagonalization method covered in class, using 
   * the eigen-decomposition functionalities provided by EJML.
   */
  @Tutorial(showSource = false, showLink = true)
  public double [][] marginalTransitionProbability(double branchLength)
  {
    double [][] result =  EJMLUtils.copyMatrixToArray(EJMLUtils.matrixExponential(eigenDecomp, branchLength));
    NumericalUtils.checkIsTransitionMatrix(result);
    for (int row = 0; row < result.length; row++)
      Multinomial.normalize(result[row]);
    return result;
  }


  public double [] stationaryDistribution()
  {
    return stationaryDistribution;
  }

  @Override
  public double[][] getRateMatrix()
  {
    return rates;
  }
  
}
