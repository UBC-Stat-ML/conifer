package conifer.ctmc.cnv;

import java.util.ArrayList;
import java.util.Arrays;
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

  
  // TODO: what are these probabilities? (In the interface it says: prior probability of each category. Let's assume we have one and return 1).
  @Override
  public List<Double> getLogPriorProbabilities()
  {
    return Arrays.asList(1.0);
  }

}
