package conifer.ctmc.cnv;

import java.util.List;

import blang.annotations.FactorComponent;
import conifer.ctmc.CTMCParameters;
import conifer.models.RateMatrixMixture;

public class CopyNumberMixture implements RateMatrixMixture
{

  public CopyNumberMixture(CopyNumberMatrix parameters)
  {
    this.parameters = parameters; 
  }
  
  @FactorComponent
  public final CopyNumberMatrix parameters;
  
  /**
   * this is not a multicat model at moment 
   */
  @Override
  public CTMCParameters getRateMatrix(int index)
  {
    return parameters; 
  }

  
   // what are these probabilities? 
  
  @Override
  public List<Double> getLogPriorProbabilities()
  {
    // TODO Auto-generated method stub
    return null;
  }

}
