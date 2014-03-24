package conifer.models;

import java.util.List;

import conifer.ctmc.CTMCParameters;

/**
 * A parameterization for a list of rate matrices with their probabilities.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public interface RateMatrixMixture
{
  public CTMCParameters getRateMatrix(int index);
  
  /**
   * @return Prior probability of each of the category.
   */
  public List<Double> getLogPriorProbabilities();

}