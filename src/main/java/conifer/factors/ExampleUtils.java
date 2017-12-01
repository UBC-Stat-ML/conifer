package conifer.factors;

import blang.core.RealConstant;
import conifer.ctmc.SimpleRateMatrix;
import conifer.models.DiscreteGammaMixture;
import conifer.models.MultiCategorySubstitutionModel;

public class ExampleUtils 
{
  
  @SuppressWarnings({ "rawtypes", "unchecked" })
  public static MultiCategorySubstitutionModel kimura(int nSites)
  { 
    return new MultiCategorySubstitutionModel(
      new DiscreteGammaMixture(
        new RealConstant(0.0), 
        new RealConstant(1.0), 
        SimpleRateMatrix.fromResource("/conifer/ctmc/kimura1980.txt"), 
        1
      ), 
      nSites
    );
  }
}
