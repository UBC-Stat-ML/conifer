package conifer.ctmc;

import blang.annotations.FactorArgument;
import blang.variables.RealVariable;



public class NoisyEmissionModel implements RateMatrixToEmissionModel
{
  private final int [] latent2observed;
  private final int nObserved;
  
  @FactorArgument
  public final RealVariable errorProbability;
  // TODO: will definitely need bounds here on the real variable
  //       if you want to allow resampling the errorPr 

  @Override
  public double[][] getMatrixStatesToObservationProbabilities()
  {
    // TODO: cache ? 
    // decide and check on a policy regarding storing transient stuff in Factors
    
    final int nLatents  = latent2observed.length;
    
    final double[][] result = new double[nLatents][nObserved];
    
    final double individualErrorPr = errorProbability.getValue()/(((double) nObserved) - 1.0);
    final double correctPr = 1.0 - errorProbability.getValue();
    
    for (int latentIndex = 0; latentIndex < nLatents; latentIndex++)
      for (int obsIndex = 0; obsIndex < nObserved; obsIndex++)
        result[latentIndex][obsIndex] = latent2observed[latentIndex] == obsIndex ? correctPr : individualErrorPr;
    
    return result;
  }

  public NoisyEmissionModel(int[] latent2observed, int nObserved,
      RealVariable errorProbability)
  {
    this.latent2observed = latent2observed;
    this.nObserved = nObserved;
    this.errorProbability = errorProbability;
  }
  
}
