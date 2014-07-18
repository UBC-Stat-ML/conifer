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


import java.io.*;
import java.util.*;

import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.junit.Assert;
import org.junit.Test;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;


// Added by Tingting
import org.jblas.Eigen;
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
  private static final double THRESHOLD = 1e-6;
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
     this.eigenDecomp = SimpleEigenDecomposition(rates); 
     this.stationaryDistribution = computeStationary();
     this.rates = rates;
     }
 
 public static SimpleEigenDecomposition SimpleEigenDecomposition(double [][] rates)
 {
   Matrix ratesMatrix = new Matrix(rates);
   EigenvalueDecomposition ed = new EigenvalueDecomposition(ratesMatrix);
   Matrix V = ed.getV();
   Matrix D =ed.getD();
   Matrix Vinverse = V.inverse();
   Matrix resultMatrix = V.times(D).times(V.inverse());
   //check if result and rates are close enough
   SimpleMatrix trueMatrix = new SimpleMatrix(rates);
   SimpleMatrix calculatedMatrix = new SimpleMatrix(resultMatrix.getArray()) ;
   if(EJMLUtils.isClose(trueMatrix, calculatedMatrix, THRESHOLD))
   {
     return new  SimpleEigenDecomposition(V, D, Vinverse);
   }else{
     throw new RuntimeException();
     }
   }
  
 
 
 
 public static class SimpleEigenDecomposition
 {
   public final Matrix V, D, Vinverse, calculatedRateMtx;

   public SimpleEigenDecomposition(Matrix v, Matrix d,
       Matrix vinverse)
   {
     V = v;
     D = d;
     Vinverse = vinverse;
     calculatedRateMtx=V.times(D).times(V.inverse());
   }
   
 }

  
  
  
  /**
   * Compute and cache the stationary (called by the constructor)
   * @return
   */
  private double[] computeStationary()
  {
    Matrix marginal = new Matrix(marginalTransitionProbability(1.0));
    Matrix exponentiatedTranspose = marginal.transpose();
   
    double [] result = null;
    final int dim = marginal.getColumnDimension();
    EigenvalueDecomposition eigenDecomp = new EigenvalueDecomposition(exponentiatedTranspose);
    // Find an eigen value equal to one (up to numerical precision)
    // result will hold the corresponding eigenvalue
    double []  realEigenvalues=eigenDecomp.getRealEigenvalues();
    double []  imageEigenvalues=eigenDecomp.getImagEigenvalues();
    double [][] eigenVectorTranspose = eigenDecomp.getV().transpose().getArray();
    for (int i = 0; i < dim; i++)
      if (NumericalUtils.isClose(0.0, imageEigenvalues[i], NumericalUtils.THRESHOLD)&& 
          NumericalUtils.isClose(1.0, realEigenvalues[i], NumericalUtils.THRESHOLD))
      {
        result = eigenVectorTranspose[i];
      }
    if (result == null)
      throw new RuntimeException("Could not find an eigenvalue equal to one. Not a proper rate matrix?");
    double sum = NumericalUtils.getNormalization(result);
    if (sum < 0.0)
    {
      for (int i = 0; i < result.length; i++)
        {
        result[i] *= -1;
        }
      }
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
    DoubleMatrix rateMtx = new DoubleMatrix(eigenDecomp.calculatedRateMtx.getArray());
    rateMtx.muli(branchLength);
    double [][] result = MatrixFunctions.expm(rateMtx).toArray2();
    NumericalUtils.checkIsTransitionMatrix(result);
   
    for (int row = 0; row < result.length; row++)
      {
        Multinomial.normalize(result[row]);
      }
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
