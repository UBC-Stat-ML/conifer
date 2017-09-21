package conifer.RandomUtils;

/**
 * Created by crystal on 2017-09-21.
 */

import blang.types.RealScalar;

import java.util.Random;

/**
 * Exponential densities.
 *
 * P is the type of parameterization.
 *
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 */
public class Exponential {

    public static double sampleUnitRateExponential(Random random)

    {

        return - Math.log(random.nextDouble());

    }


    public static double sampleExponential(Random random, double rate)

    {

        return sampleUnitRateExponential(random) / rate;

    }

    public static double generate(Random random, double rate){
        return sampleExponential(random, rate);
    }

  

}