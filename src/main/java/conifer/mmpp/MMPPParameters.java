package conifer.mmpp;

import conifer.ctmc.CTMC;
import conifer.ctmc.CTMCParameters;
import conifer.ctmc.RateMatrixToEmissionModel;

/**
 * Created by crystal on 2016-04-19.
 */
public interface MMPPParameters extends CTMCParameters{

    public double[] getIntensities();
}


