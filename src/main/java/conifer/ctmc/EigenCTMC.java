package conifer.ctmc;

import org.ejml.data.DenseMatrix64F;
import org.ejml.data.Eigenpair;
import org.ejml.factory.DecompositionFactory;
import org.ejml.factory.EigenDecomposition;
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
    
    double [] result = null;
    final int dim = marginal.numCols();
    EigenDecomposition<DenseMatrix64F> eigenDecomp = DecompositionFactory.eig(exponentiatedTranspose.numCols(), true);
    eigenDecomp.decompose(exponentiatedTranspose.getMatrix());
    // Find an eigen value equal to one (up to numerical precision)
    // result will hold the corresponding eigenvalue
    for (int i = 0; i < dim; i++)
      if (eigenDecomp.getEigenvalue(i).isReal() && 
          NumericalUtils.isClose(1.0, eigenDecomp.getEigenvalue(i).getReal(), NumericalUtils.THRESHOLD))
        result = eigenDecomp.getEigenVector(i).data;
    if (result == null)
      throw new RuntimeException("Could not find an eigenvalue equal to one. Not a proper rate matrix?");
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
    
    // For small branch lengths, some small negative values can appear. Set those to a small positive value instead.
    for (int i = 0; i < result.length; i++)
      for (int j = 0; j < result.length; j++)
        if (i != j && result[i][j] < 0.0)
        {
          if (result[i][j] + NumericalUtils.THRESHOLD < 0.0)
            throw new RuntimeException();
          result[i][j] = -result[i][j];
        }
    
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
