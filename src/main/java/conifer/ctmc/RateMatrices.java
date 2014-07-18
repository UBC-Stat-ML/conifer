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

  public static SimpleRateMatrix rateMtxModel(String selectedRateMtx)
  {
    SimpleRateMatrix rateMtx;
    RateMatrices result = new RateMatrices();
    rateMtx=null;
    
    if (selectedRateMtx == null) {
      throw new IllegalArgumentException("selectedRateMtx is null!");
  }

  final RateMtxNames something = RateMtxNames.fromString(selectedRateMtx);

  if (something == null) {
     return result.kimura1980();
  }

  switch(something) {
      case KIMURA1980:
          return result.kimura1980();
      case ACCORDANCE:
          return result.accordance();
      case PAIR:
         // return result.pair()
      case POLARITY:
          return result.polarity();
      case POLARITYSIZE:
          return result.polaritySize();
      
      default:
          return result.kimura1980();  
  }
    
//    {
//      Random rand = new Random();
//      int size =400;
//      rateMtx=result.randomGTR(rand, size);
//    }
     
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
