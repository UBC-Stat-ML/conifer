package conifer.ctmc;

import java.util.Random;

import conifer.ctmc.expfam.RateMtxNames;
import conifer.ctmc.expfam.SerializedExpFamMixture;

import bayonet.distributions.Multinomial;



public class RateMatrices
{
  public static SimpleRateMatrix kimura1980()
  {
    return SimpleRateMatrix.fromResource("/conifer/ctmc/kimura1980.txt");
  }
  
  public static SimpleRateMatrix accordance()
  {
    return SimpleRateMatrix.fromResource("/conifer/ctmc/accordance.txt");
  }
  
  public static SimpleRateMatrix polarity()
  {
    return SimpleRateMatrix.fromResource("/conifer/ctmc/polarity.txt");
  }

  public static SimpleRateMatrix polaritySize()
  {
    return SimpleRateMatrix.fromResource("/conifer/ctmc/polaritySize.txt");
  }

  public static SimpleRateMatrix polaritySizeGTR()
  {
    return SimpleRateMatrix.fromResource("/conifer/ctmc/polaritySizeGTR.txt");
  }
  
  public static SimpleRateMatrix rateMtxModel(final RateMtxNames selectedRateMtx)
  {
    SimpleRateMatrix rateMtx;
    RateMatrices result = new RateMatrices();
    rateMtx=null;
    
    if (selectedRateMtx == null) {
      throw new IllegalArgumentException("selectedRateMtx is null!");
  }
    return selectedRateMtx.getRateMtx();     
  }
  
  /**
   * Generate uniform pi and thetas, normalize the pi, and create a 
   * GTR from this. 
   * Note: not so natural of a distribution, but simple to implement 
   * and useful for testing.
   * Warning: not normalized to have an expected number of change to be one
   * for branch length one.
   * 
   * Also destructively writes the true statio in the stationary argument.
   * 
   * @param rand
   * @return
   */
  public static SimpleRateMatrix randomGTR(Random rand, int size, double [] stationary)
  {
    if (stationary == null)
      stationary = new double[size];
    double [] thetas = new double[size*(size-1)/2];
    for (int i = 0; i < stationary.length; i++)
      stationary[i] = rand.nextDouble();
    Multinomial.normalize(stationary);
    for (int i = 0; i < thetas.length; i++)
      thetas[i] = rand.nextDouble();
    double [][] matrix = RateMatrixUtils.gtrFromOverParam(stationary, thetas, size);
    return new SimpleRateMatrix(matrix, null);
  }
  
  /**
   * Generate uniform pi and thetas, normalize the pi, and create a 
   * GTR from this. 
   * Note: not so natural of a distribution, but simple to implement 
   * and useful for testing.
   * Warning: not normalized to have an expected number of change to be one
   * for branch length one.
   * 
   * @param rand
   * @return
   */
  public static SimpleRateMatrix randomGTR(Random rand, int size)
  {
    return randomGTR(rand, size, null);
  }
}
