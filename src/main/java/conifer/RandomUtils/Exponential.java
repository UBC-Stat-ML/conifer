package conifer.RandomUtils;

/**
 * Created by crystal on 2017-09-21.
 */

import blang.core.RealVar;
import blang.types.RealScalar;

import java.util.Random;

/**
 * Exponential densities.
 *
 * P is the type of parameterization.
 *
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 */
public class Exponential<P> {

    public final blang.core.RealVar realization;
    /**
     * The parameter of this exponential density.
     */

    public final P parameters;


    public Exponential(RealVar realization, P rateParameterization) {

        this.realization = realization;
        this.parameters = rateParameterization;
    }

    public static class RateParameterization
    {
      /**
       * 
       */
      
      public final blang.core.RealVar rate;
      
      /**
       * 
       * @param rate
       */
      public RateParameterization(blang.core.RealVar rate)
      {
        this.rate = rate;
      }

      /**
       * @return The rate.
       */
    
      public double getRate()
      {
        return rate.doubleValue();
      }
    
    }

    public static Exponential<RateParameterization> on(blang.core.RealVar realization)
    {
        return new Exponential<RateParameterization>(realization, new RateParameterization(new RealScalar(1.0)));
    }

    public static Exponential<RateParameterization> on(RealScalar realization){
        return new Exponential<RateParameterization>(realization, new RateParameterization(new RealScalar(1.0)));
    }


    public static interface Parameters
    {
        /**
         * @return The parameter, transformed back into a rate.
         */
        public double getRate();
    }


    public static double sampleUnitRateExponential(Random random)

    {

        return - Math.log(random.nextDouble());

    }


    public static double sampleExponential(Random random, double rate)

    {

        return sampleUnitRateExponential(random) / rate;

    }
    
    public static Exponential<RateParameterization> newExponential()
    {
      return Exponential.on(new RealScalar(1.0));
    }
    

    public static double generate(Random random, double rate){
        return sampleExponential(random, rate);
    }


  
}