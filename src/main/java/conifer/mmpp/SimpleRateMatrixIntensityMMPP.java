package conifer.mmpp;

import conifer.ctmc.EigenCTMC;
import conifer.ctmc.RateMatrixToEmissionModel;
import conifer.ctmc.SimpleRateMatrix;

/**
 * Created by crystal on 2016-04-20.
 */
public class SimpleRateMatrixIntensityMMPP extends SimpleRateMatrix implements MMPPParameters
{

    private final double [] intensities;

    public SimpleRateMatrixIntensityMMPP(double[][] rateMatrix, RateMatrixToEmissionModel emissionModel, double[] intensities){
        super(rateMatrix, emissionModel);
        this.intensities=intensities;
    }

    @Override
    public double[] getIntensities() {
        return this.intensities;
    }
}
