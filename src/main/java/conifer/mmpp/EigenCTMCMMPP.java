package conifer.mmpp;

import bayonet.distributions.Multinomial;
import bayonet.math.NumericalUtils;
import conifer.ctmc.EigenCTMC;
import conifer.ctmc.RateMatrixUtils;

/**
 * Created by crystal on 2016-04-20.
 */
public class EigenCTMCMMPP extends EigenCTMC implements MMPP {

    /**
     * Note: if the RateMatrix is changed in place,
     * these changes will not be mirrored by this class.
     * <p>
     * It should be recreated each time a likelihood
     * calculation is performed.
     *
     * @param rates
     */
    private final double[] intensities;
    public EigenCTMCMMPP(double[][] rates, double[] intensities) {
        super(rates);
        this.intensities=intensities;

    }

    @Override
    public double [][] marginalTransitionProbability(double [][] rateMatrix, double [] intensities, double T){

        double [][] rateMatrixMinusLamda = RateMatrixUtils.rateMatrixMinusDiagIntensity(rateMatrix, intensities);
        double [][] result = new double [rateMatrixMinusLamda[0].length][rateMatrixMinusLamda[0].length];

        result = RateMatrixUtils.marginalTransitionMtx(rateMatrixMinusLamda, T, false);

        //NumericalUtils.checkIsTransitionMatrix(result);

        //for (int row = 0; row < result.length; row++)
        //{
        //    Multinomial.normalize(result[row]);
        //}
        return result;

    }


    @Override
    public double[] getIntensities() {
        return this.intensities;
    }

    @Override
    public double[][] getVirtualRateMatrix() {
        return RateMatrixUtils.fillVirtualRateMatrixForMMPPs(super.getRateMatrix(), intensities);
    }
}
