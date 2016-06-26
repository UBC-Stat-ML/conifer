package conifer.ctmc;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import bayonet.distributions.Multinomial;
import bayonet.math.NumericalUtils;
import conifer.mmpp.MMPP;
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
   * Fill in place the last element for the virtual state of MMPPs.
   * The last element in each row is the Poisson intensities for each state
   * The last row is filled up with zeros since the virtual state for Poisson event since we assume it is an absorbing state
   * @author crystal
   * @param intensities, rate
   */
   public static double [][] fillVirtualRateMatrixForMMPPs(final double [][] rates, double [] intensities)
   {
       int nStates=intensities.length;
       int nStatesIncludeVirtualState = nStates+1;
       if(intensities.length!=rates[0].length){
           throw new RuntimeException("The number of States in the rate matrix should be the same as the dimension of the Poisson intensities");
       }
       double [][] virtualRateMtx = new double[nStatesIncludeVirtualState][nStatesIncludeVirtualState];
       for(int i=0; i< nStates; i++)
       {
           for(int j=0;j<nStates; j++)
           {
               virtualRateMtx[i][j]=rates[i][j];
           }
       }

       for(int j=0;j<nStates;j++){
           virtualRateMtx[j][nStates]=intensities[j];
           virtualRateMtx[nStates][j]=0;

       }
       virtualRateMtx[nStates][nStates]=0;

       fillRateMatrixDiagonalEntriesAllowZero(virtualRateMtx);
       return virtualRateMtx;
   }

    /**
     * Fill in place the diagonals to be equal for each row to the negative of the
     * off diagonals.
     * The reason why we allow zero here is that in the virtual rate matrix for MMPPs,
     * the last row is filled with zeros for the virtual state.
     * @param rate
     */
    public static void fillRateMatrixDiagonalEntriesAllowZero(final double [][] rate)
    {
        int size = rate.length;
        for (int i = 0; i < size; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < size; j++)
                if (i!=j)
                    sum += rate[i][j];
            rate[i][i] = -sum;
        }
    }

    // added by Crystal to test correctness of fillVirtualRateMatrixForMMPPs
    public static void main(String [] args){
        double [] intensities = {1,2,3};
        double [][] rates= { {0, 1, 2} , { 4,0, 5}, {3,2, 0} };
        double [][] result = fillVirtualRateMatrixForMMPPs(rates, intensities);
        System.out.println(Arrays.deepToString(result));
        double [][] embeddedMMPPs = getJumpProcessMMPP(result);
        System.out.println(Arrays.deepToString(embeddedMMPPs));
    }

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

    public static double [][] intensityToDiagMatrix(double [] intensities)
    {
        //create a diagonal array with intensities as diagonal elements
        double [][] lambda= new double[intensities.length][intensities.length];
        for(int i=0;i<intensities.length;i++)
            lambda[i][i] = intensities[i];
        return lambda;

    }

    public static double [][] rateMatrixMinusDiagIntensity(double [][] rateMatrix, double [] intensities)
    {
        double [][] lambda = intensityToDiagMatrix(intensities);
        DoubleMatrix lambdaMatrix = new DoubleMatrix(lambda);
        int nStates = rateMatrix[0].length;
        //Then we obtain Q-lambda for MMPP
        DoubleMatrix QStar = new DoubleMatrix(rateMatrix);
        QStar = QStar.sub(lambdaMatrix);
        return QStar.toArray2();
    }

    public static double [][] rateMatrixMinusDiagIntensity(MMPP mmpp){
        double [][] rateMatrix = mmpp.getRateMatrix();
        double [] intensities = mmpp.getIntensities();
        double [][] QStar= rateMatrixMinusDiagIntensity(rateMatrix, intensities);
        return QStar;
    }

    /**
     * if the CTMC has n state, the virtual rate matrix for MMPP is (n+1)*(n+1)
     * the last row of this matrix is filled with all zeros
     * the getJumpProcessMMPP will extract off-diagonal elements of the first n rows
     * the parameter rate here is the virtual rate matrix
     * @param virtualRate
     * @return
     */

  public static double [][] getJumpProcessMMPP(double [][] virtualRate)
  {
      int nrow = (virtualRate[0].length-1);
      int ncolumn = virtualRate[0].length;
      double [][] result = new double[nrow][ncolumn];

      for(int i=0; i<nrow;i++){
          for(int j=0; j<ncolumn;j++)
          {
              if(i!=j)
                  result[i][j]=virtualRate[i][j];

          }
          Multinomial.normalize(result[i]);
      }
      return result;
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


    //ToDo: this should go to a MatrixUtils in bayonet
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
