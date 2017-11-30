package conifer.RandomUtils;

import java.util.Random;

import bayonet.math.SpecialFunctions;
import blang.core.RealVar;
import blang.types.RealScalar;
import conifer.RandomUtils.Normal.Parameters;

/**
 * Created by crystal on 2017-09-21.
 * @param <P>
 */
public class Normal<P> implements PriorOfWeights{

    public static double quantile(double z, double m, double sd)
    {
        return m + Math.sqrt(2.0) * sd * SpecialFunctions.inverseErf(2.0 * z - 1.0);
    }
    
    
    public RealVar realization;
    
    public final P parameters;
    
    public static interface Parameters
    {
      public double getMean();
      public double getVariance();
    }
    
    public Normal(RealVar realization, P parameters)
    {
      this.realization = realization;
      this.parameters = parameters;
    }
    
    public Normal(P parameters){
    	this.parameters = parameters;
    }
    
    public final static double LOG_INV_SQRT_2_PI = -Math.log(Math.sqrt(2*Math.PI));
    
    /**
     * The product of  normal densities for n iid observations x_i
     * @param mean
     * @param var
     * @param sum \sum_i x 
     * @param sumSq \sum_i x^2
     * @param n
     * @return
     */
    public static double logProb(double mean, double var, double sum, double sumSq, double n) 
    {
      return -0.5*(sumSq - 2*mean*sum + mean*mean)/var + n*(LOG_INV_SQRT_2_PI - 0.5*Math.log(var));
    }
    
    public static double logDensity(double x, double mean, double var)
    {
      return logProb(mean, var, x, x*x, 1);
    }
    
    /**
     * Box-Muller Transformation
     * 
     * @param random
     * @return N(0,1) sample
     */
    public static double generateStandardNormal(Random random) 
    {
      double x1 = random.nextDouble(), x2 = random.nextDouble();
      double z = Math.sqrt(-2*Math.log(x1))*Math.cos(2*Math.PI*x2);
      return z;
    }
    
    public static double generate(Random random, double mean, double var) 
    {
      return generateStandardNormal(random) * Math.sqrt(var) + mean;
    }
    
    
 
    public double logDensity()
    {
      return logDensity(realization.doubleValue(), ((Parameters) parameters).getMean(), ((Parameters) parameters).getVariance());
    }
    
    @Override
    public double logDensity(double x)
    {
      return logDensity(x, ((Parameters) parameters).getMean(), ((Parameters) parameters).getVariance());
    }
    

   
    public RealVar getRealization()
    {
      return realization;
    }

  
   // public void generate(Random random)
   // {
   //   realization.setValue(generate(random, ((Parameters) parameters).getMean(), ((Parameters) parameters).getVariance()));
   // }
    
    public static class MeanVarianceParameterization implements Parameters
    {
      
      public final RealVar mean;
      

      public final RealVar variance;
      
      public MeanVarianceParameterization(RealVar mean, RealVar var)
      {
        this.mean = mean;
        this.variance = var;
      }

      @Override     
      public double getMean()
      {
        return mean.doubleValue();
      }

      @Override
      public double getVariance()
      {
        return variance.doubleValue();
      }
    }
    
    public static class MeanPrecisionParameterization implements Parameters
    {
      
      public final RealVar mean; 
      
   
      public final RealVar precision;
        
      public MeanPrecisionParameterization(RealVar mean, RealVar precision)
      {
          this.mean = mean; 
          this.precision = precision;
      }
      
      @Override
      public double getMean()
      {
          return mean.doubleValue();
      }

      @Override
      public double getVariance()
      {
          return (1 / precision.doubleValue());
      }
        
    }
    
    /* Syntactic sugar/method chaining */
    
    public static Normal<MeanVarianceParameterization> on(RealVar realization)
    {
      return new Normal<MeanVarianceParameterization>(realization, new MeanVarianceParameterization(new RealScalar(0.0), new RealScalar(1.0)));
    }
    
    public static Normal<MeanVarianceParameterization> newNormal()
    {
      return Normal.on(new RealScalar(0.0));
    }
      
    public Normal<MeanPrecisionParameterization> withMeanPrecision(RealVar mean, RealVar precision)
    {
        return new Normal<MeanPrecisionParameterization>(realization, new MeanPrecisionParameterization(mean, precision)); 
    }

	

	@Override
	public double gradientOflogDensity(double x) {
		double result = -(x- ((Parameters) parameters).getMean())/((Parameters) parameters).getVariance();
		return result;
	}

	@Override
	public double generate(Random random, Object parameterization) {
		return generate(random, ((Parameters) parameters).getMean(), ((Parameters) parameters).getVariance());
	}

}
