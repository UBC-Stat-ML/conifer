package conifer.models;

import bayonet.distributions.Multinomial;



public class RateMatrixMixtureUtils
{
  public static int nCategories(RateMatrixMixture mixture)
  {
    return mixture.getLogPriorProbabilities().size();
  }

//  public static double[] priorProbabilities(RateMatrixMixture mixture)
//  {
//    double[] result = new double[nCategories(mixture)];
//    for (int i = 0; i < nCategories(mixture); i++)
//      result[i] = mixture.getLogPriorProbabilities().get(i);
//    Multinomial.expNormalize(result);
//    return result;
//  }
}
