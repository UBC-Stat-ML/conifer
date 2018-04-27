package conifer.models;

import java.util.List;


import bayonet.distributions.Multinomial;
import blang.core.RealVar;

import com.google.common.collect.Lists;

import conifer.Utils;
import conifer.ctmc.CTMC;
import conifer.ctmc.CTMCParameters;
import conifer.ctmc.InvariantCTMCParameters;
import conifer.ctmc.SimpleRateMatrix;



public class DiscreteGammaMixture implements RateMatrixMixture
{

  public final RealVar invariantSiteProbability;

  public final RealVar shapeParameter;

  public final CTMCParameters baseRateMatrix;
  
  private final int nPositiveCategories;
  
  public DiscreteGammaMixture(
      RealVar invariantSiteProbability,
      RealVar shapeParameter,
      CTMCParameters baseRateMatrix,
      int nPositiveCategories)
  {
    this.invariantSiteProbability = invariantSiteProbability;
    this.shapeParameter = shapeParameter;
    this.baseRateMatrix = baseRateMatrix;
    this.nPositiveCategories = nPositiveCategories;
  }
  
  @Override
  public CTMCParameters getRateMatrix(int index)
  {
    double rate = computeDiscreteGammaRates(invariantSiteProbability.doubleValue(), nPositiveCategories, shapeParameter.doubleValue()).get(index);
    if (rate == 0.0)
    {
      CTMC base = baseRateMatrix.getProcess();
      return new InvariantCTMCParameters(base.stationaryDistribution(), baseRateMatrix.getEmissionModel());
    }
    else
    {
      double [][] baseMatrix = baseRateMatrix.getRateMatrix();
      final int matrixSize = baseMatrix.length;
      double [][] scaled = new double[matrixSize][matrixSize];
      if (rate != 0.0)
        for (int i = 0; i < matrixSize; i++)
          for (int j = 0; j < matrixSize; j++)
            scaled[i][j] = baseMatrix[i][j] * rate;
      return new SimpleRateMatrix(scaled, baseRateMatrix.getEmissionModel());
    }
  }

  @Override
  public List<Double> getLogPriorProbabilities()
  {
    return computeDiscreteGammaLogPriorProbabilities(invariantSiteProbability.doubleValue(), nPositiveCategories, shapeParameter.doubleValue());
  }
  
  public static List<Double> computeDiscreteGammaRates(double invariantCategoryPr, int nPositiveCategories, double shapeParameter)
  {
    if (shapeParameter <= 0.0 || invariantCategoryPr < 0.0 || invariantCategoryPr > 1.0)
      throw new RuntimeException();
    final List<Double> result = Lists.newArrayList();
    boolean useInvar = invariantCategoryPr > 0.0;
    
    double remainingMass = 1.0;
    if (useInvar)
    {
      remainingMass =  1.0 - invariantCategoryPr;
      result.add(0.0);
    }
    
    double [] rates = new double[nPositiveCategories];
    for (int i = 0; i < nPositiveCategories; i++)
      rates[i] = Utils.gammaQuantile((2.0 * i + 1.0) / (2.0 * nPositiveCategories), shapeParameter, 1.0 / shapeParameter); // from BEAST
    final double mean = remainingMass * Multinomial.getNormalization(rates) / nPositiveCategories;
    
    for (int i = 0; i < nPositiveCategories; i++)
    {
      double normalizedRate = rates[i] / mean;
      result.add(normalizedRate);
    }
    
    return result;
  }
  
  public static List<Double> computeDiscreteGammaLogPriorProbabilities(double invariantCategoryPr, int nPositiveCategories, double shapeParameter)
  {
    if (shapeParameter <= 0.0 || invariantCategoryPr < 0.0 || invariantCategoryPr > 1.0)
      throw new RuntimeException();
    final List<Double> result = Lists.newArrayList();
    boolean useInvar = invariantCategoryPr > 0.0;
    
    double remainingMass = 1.0;
    if (useInvar)
    {
      remainingMass =  1.0 - invariantCategoryPr;
      result.add(Math.log(invariantCategoryPr));
    }
    
    for (int i = 0; i < nPositiveCategories; i++)
      result.add(Math.log(remainingMass/nPositiveCategories));
    
    return result;
  }
  
  public static void main(String [] args)
  {
    for (double i = 0.1; i < 10; i *= 2)
      System.out.println(computeDiscreteGammaRates(0.5, 3, i));
  }

}
