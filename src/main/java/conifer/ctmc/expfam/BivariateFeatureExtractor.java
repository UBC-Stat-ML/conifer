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
  /**
   * Note: the implementor does not to have to worry about maintaining
   * an unordered representation of the two states: for any two distinct 
   * states {s1, s2}, features will be 
   * extracted for exactly one of (s1, s2) or (s2, s1). The features for 
   * the permuted one will be set automatically to be the same set of bivariate
   * features to maintain reversibility.
   * 
   * @param counts
   * @param state1
   * @param state2
   */
  @SuppressWarnings("rawtypes")
  public void extract(Counter counts, S state1, S state2);
}