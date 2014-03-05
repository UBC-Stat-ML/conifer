package conifer.ctmc;



public class JukeCantorRateMatrix implements RateMatrix
{
  private final double [][] matrix;
  
  public JukeCantorRateMatrix(int nCharacters)
  {
    if (nCharacters <= 1)
      throw new RuntimeException();
    matrix = new double[nCharacters][nCharacters];
    double entry = 1.0 / (nCharacters - 1);
    for (int i = 0; i < nCharacters; i++)
      for (int j = 0; j < nCharacters; j++)
        if (i != j)
          matrix[i][j] = entry;
    RateMatrixUtils.fillRateMatrixDiagonalEntries(matrix);
  }

  @Override
  public double[][] getMatrix()
  {
    return matrix;
  }

}
