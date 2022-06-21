package conifer.ctmc;

import blang.types.StaticUtils;

public class NoisyEmissionModel implements RateMatrixToEmissionModel
{
  private final int [] latent2observed;
  private final int nObserved;

  public final blang.core.RealVar errorProbability;

  @Override
  public double[][] getMatrixStatesToObservationProbabilities(double observationAnnealing)
  {
    if (observationAnnealing < 1.0 && errorProbability.doubleValue() == 0.0)
      throw new RuntimeException();
    
    if (errorProbability.doubleValue() < 0.0 || errorProbability.doubleValue() > 1.0)
      StaticUtils.invalidParameter();
        
    final int nLatents  = latent2observed.length;
    final double[][] result = new double[nLatents][nObserved];
    
    double individualErrorPr = errorProbability.doubleValue()/(((double) nObserved) - 1.0);
    double correctPr = 1.0 - errorProbability.doubleValue();
    
    if (observationAnnealing < 1.0) 
    {
      final double individualErrorPrAnnealedUnnorm = Math.pow(individualErrorPr, observationAnnealing);
      final double correctPrAnnealedUnnorm = Math.pow(correctPr, observationAnnealing);
      
      final double norm = correctPrAnnealedUnnorm + (nObserved - 1.0) * individualErrorPrAnnealedUnnorm;
      
      individualErrorPr = individualErrorPrAnnealedUnnorm / norm;
      correctPr = correctPrAnnealedUnnorm / norm;
    }
      
    for (int latentIndex = 0; latentIndex < nLatents; latentIndex++)
      for (int obsIndex = 0; obsIndex < nObserved; obsIndex++)
        result[latentIndex][obsIndex] = latent2observed[latentIndex] == obsIndex ? correctPr : individualErrorPr;
    
    return result;
  }

  public NoisyEmissionModel(int[] latent2observed, int nObserved,
      blang.core.RealVar errorProbability)
  {
    this.latent2observed = latent2observed;
    this.nObserved = nObserved;
    this.errorProbability = errorProbability;
  }
  
}
