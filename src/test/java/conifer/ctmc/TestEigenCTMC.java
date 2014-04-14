package conifer.ctmc;

import java.util.Arrays;
import java.util.Random;


import org.ejml.simple.SimpleMatrix;
import org.junit.Test;

import bayonet.math.NumericalUtils;



public class TestEigenCTMC
{
  @Test
  public void testStationary()
  {
    Random rand = new Random(1);
    for (int j = 0; j < 10000; j++)
    {
      int size = 10;
      double [] trueStatio = new double[size];
      SimpleRateMatrix rateM = RateMatrices.randomGTR(rand, size, trueStatio);
      EigenCTMC ctmc = new EigenCTMC(rateM.getRateMatrix());
      double [] statio = new SimpleMatrix(size,1,true,ctmc.stationaryDistribution()).getMatrix().data;
      for (int i = 0; i < size; i++)
        NumericalUtils.checkIsClose(statio[i], trueStatio[i]);
    }
  }
}
