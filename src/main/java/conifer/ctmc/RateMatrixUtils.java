package conifer.ctmc;

import bayonet.math.NumericalUtils;


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

  /**
   * Check off diagonals are non-negative and  sum to the negative of diagonals.
   * @param rateMatrix
   */
  public static void checkValidRateMatrix(RateMatrix rateMatrix)
  {
    double [][] rates = rateMatrix.getMatrix();
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
}
