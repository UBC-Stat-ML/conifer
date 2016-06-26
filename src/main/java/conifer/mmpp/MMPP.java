package conifer.mmpp;

import conifer.ctmc.CTMC;

/**
 * Created by crystal on 2016-04-19.
 */
public interface MMPP extends CTMC
{
    public double [] getIntensities();
    public double [][] getVirtualRateMatrix();
    public double [][] marginalTransitionProbability(double [][] rateMatrix, double [] intensities, double branchLength);
}
