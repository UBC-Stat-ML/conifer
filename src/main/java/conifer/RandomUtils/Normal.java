package conifer.RandomUtils;

import bayonet.math.SpecialFunctions;

/**
 * Created by crystal on 2017-09-21.
 */
public class Normal {

    public static double quantile(double z, double m, double sd)
    {
        return m + Math.sqrt(2.0) * sd * SpecialFunctions.inverseErf(2.0 * z - 1.0);
    }

}
