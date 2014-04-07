package conifer.ctmc.expfam.features;

import briefj.collections.Counter;
import conifer.ctmc.expfam.BivariateFeatureExtractor;



public class IdentityBivariate<T> implements BivariateFeatureExtractor<T>
{

  @SuppressWarnings({ "rawtypes", "unchecked" })
  @Override
  public void extract(Counter counts, T state1, T state2)
  {
    counts.incrementCount("BivariateIdentity[" + state1 + "," + state2 + "]", 1.0);
  }
  

}
