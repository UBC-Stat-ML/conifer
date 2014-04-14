package conifer.ctmc.expfam.features;

import briefj.collections.Counter;
import conifer.ctmc.expfam.UnivariateFeatureExtractor;



public class IdentityUnivariate<T> implements UnivariateFeatureExtractor<T>
{

  @SuppressWarnings({ "rawtypes", "unchecked" })
  @Override
  public void extract(Counter counts, T state)
  {
    counts.incrementCount("UnivariateIdentity[" + state.toString() + "]", 1.0);
  }

}
