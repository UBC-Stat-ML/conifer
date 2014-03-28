package conifer.ctmc.expfam;

import briefj.collections.Counter;

/**
 * Features concerning a pair of states
 * @author bouchard
 *
 * @param <S>
 */
public interface BivariateFeatureExtractor<S>
{
  public void extract(Counter counts, S state1, S state2);
}