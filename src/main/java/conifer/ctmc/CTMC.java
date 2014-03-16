package conifer.ctmc;

import org.ejml.data.Eigenpair;
import org.ejml.ops.EigenOps;
import org.ejml.simple.SimpleMatrix;

import tutorialj.Tutorial;
import bayonet.distributions.Multinomial;
import bayonet.math.EJMLUtils;
import bayonet.math.EJMLUtils.SimpleEigenDecomposition;
import bayonet.math.NumericalUtils;

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
public class CTMC
{
  private final SimpleEigenDecomposition eigenDecomp;
  private final double [] stationaryDistribution;

  /**
   * Note: if the RateMatrix is changed in place,
   * these changes will not be mirrored by this class.
   * 
   * It should be recreated each time a likelihood 
   * calculation is performed.
   * @param rateMatrix
   */
  public CTMC(RateMatrix rateMatrix)
  {
    RateMatrixUtils.checkValidRateMatrix(rateMatrix);
    this.eigenDecomp = EJMLUtils.simpleEigenDecomposition(new SimpleMatrix(rateMatrix.getMatrix()));
    this.stationaryDistribution = computeStationary();
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

  /**
   * Fill in ``stationaryDistribution()``.
   * 
   * Show how you can reduce the problem of finding the stationary distribution
   * of a CTMC to a eigenvalue problem. Again, solve this problem using
   * EJML.
   */
  @Tutorial(showSource = false, showLink = true)
  public double [] stationaryDistribution()
  {
    return stationaryDistribution;
  }
  
}
