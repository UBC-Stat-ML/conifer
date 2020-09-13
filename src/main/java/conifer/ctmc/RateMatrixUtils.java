package conifer.ctmc;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import bayonet.distributions.Multinomial;
import bayonet.math.NumericalUtils;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import java.util.Arrays;


/**
 * Utilities for rate matrices.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class RateMatrixUtils
{
  /**
   * Fill in place the diagonals to be equal for each row to the negative of the
   * off diagonals. 
   * @param rate
   */
  public static void fillRateMatrixDiagonalEntries(final double [][] rate)
  {
    int size = rate.length;
    for (int i = 0; i < size; i++)
    {
      double sum = 0.0;
      for (int j = 0; j < size; j++)
        if (i!=j)
          sum += rate[i][j];
      if (rate[i][i] != 0.0)
        throw new RuntimeException();
      rate[i][i] = -sum;
    }
  }
  
  public static double [][] getJumpProcess(double [][] rate)
  {
    final int size = rate.length;
    double [][] result = new double[size][size];
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
        if (i!=j)
          result[i][j] = rate[i][j];
      Multinomial.normalize(result[i]);
    }
    return result;
  }


  /**
   * Check off diagonals are non-negative and  sum to the negative of diagonals.
   */
  public static void checkValidRateMatrix(double [][] rates)
  {
    int size = rates.length;
    for (int row = 0; row < size; row++)
    {
      if (rates[row].length != size)
        throw new RuntimeException("Rate matrices need to be square matrices");
      double sum = 0;
      for (int col = 0; col < size; col++)
        if (col != row)
        {
          final double cur = rates[row][col];
          if (cur < 0)
            throw new RuntimeException();
          sum += cur;
        }
      NumericalUtils.checkIsClose(-rates[row][row], sum);
    }
  }
  
  /**
   * Build a GTR from parameterized vectors
   * 
   * For example, with DNA (n=4), use an array stat[] of length 4 (summing to one),
   * and an array subRates[] of length 6
   * 
   * Note: this does not take care of normalizing it to an expected number of change
   * equal to one.
   */
  public static double [][] gtrFromOverParam(double [] stat, double [] subRates, int n)
  {
    if (stat.length != n || subRates.length != n*(n-1)/2)
      throw new RuntimeException();
    
    NumericalUtils.checkIsProb(stat);
    
    double [][] result = new double[n][n];
    
    int cur = 0;
    for (int col = 0; col < n; col++)
      for (int row = 0; row < n; row++)
        if (col < row)
          result[row][col] = subRates[cur++] * stat[col];
          
    cur = 0;
    for (int row = 0; row < n; row++)
      for (int col = 0; col < n; col++)
        if (col > row)
          result[row][col] = subRates[cur++] * stat[col];
    
    
    fillRateMatrixDiagonalEntries(result);
    
    return result;
  }

    // The following code is added by zhaott0416@gmail.com
    // In order to use diagonalization for the matrix exponential algorithm
    public static enum MatrixExponentialAlgorithm
    {
        DIAGONALIZATION {
            @Override
            public double[][] marginalTransitionMtx(double[][] rate, double t)
            {
                if (Math.abs(t) < 1e-10) // protect against numerical problems
                    return Matrix.identity(rate.length, rate.length).getArray();

                final Matrix M = new Matrix(rate).times(t);
                double [][] p =  exp(M).getArray();
                for (int i = 0; i < p.length; i++)
                {
                    NumericalUtils.checkIsProb(p[i]);
                    Multinomial.normalize(p[i]);
                }
                return p;
            }
        },
        BLAS {
            @Override
            public double[][] marginalTransitionMtx(double[][] rate, double t)
            {
                DoubleMatrix rateMtx = new DoubleMatrix(rate);
                rateMtx.muli(t);
                return MatrixFunctions.expm(rateMtx).toArray2();
            }
        };
        public abstract double [][] marginalTransitionMtx(final double [][] rate, final double t);
    }

    public static MatrixExponentialAlgorithm defaultsMatrixExponentialAlgorithm(boolean useDiag){

        MatrixExponentialAlgorithm result;
        if(useDiag)
            result = MatrixExponentialAlgorithm.DIAGONALIZATION;
        else
            result = MatrixExponentialAlgorithm.BLAS;

        return result;
    }

    public static double [][] marginalTransitionMtx(final double [][] rate, final double t){

        return MatrixExponentialAlgorithm.BLAS.marginalTransitionMtx(rate, t);
    }


    public static double [][] marginalTransitionMtx(final double [][] rate, final double t, MatrixExponentialAlgorithm method)
    {
        return method.marginalTransitionMtx(rate, t);
    }

    public static double [][] marginalTransitionMtx(final double [][] rate, final double t, boolean useDiag){
        return defaultsMatrixExponentialAlgorithm(useDiag).marginalTransitionMtx(rate, t);
    }

    public static Matrix exp(final Matrix V, final Matrix Vinv, final Matrix D)
    {
        final int size = V.getColumnDimension();
        final Matrix expD = new Matrix(size, size);
        for (int i = 0; i < size; i++)
            expD.set(i, i, Math.exp(D.get(i, i)));
        Matrix result = V.times(expD).times(Vinv);
        for (int row = 0; row < V.getRowDimension(); row++)
        {
            double sum = 0.0;
            for (int col = 0; col < V.getColumnDimension(); col++)
                sum += result.get(row, col);
            NumericalUtils.checkIsClose(sum, 1.0);
        }
        return result;
    }

    public static Matrix exp(final Matrix M)
    {
        final EigenvalueDecomposition ed = new EigenvalueDecomposition(M);
        final Matrix V = ed.getV();
        return exp(V, V.inverse(), ed.getD());
    }

}
