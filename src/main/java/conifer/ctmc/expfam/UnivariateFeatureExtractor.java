package conifer.ctmc.expfam;

import briefj.collections.Counter;

/**
 * Features concerning a single state (e.g. for stationary or perhaps total rates)
 * @author bouchard
 *
 * @param <S>
 */
public interface UnivariateFeatureExtractor<S>
{
  public void extract(Counter counts, S state);
}