package conifer.ctmc;

import java.util.Random;

import org.ejml.simple.SimpleMatrix;
import org.junit.Test;

import bayonet.math.EJMLUtils;



/**
 * Test the agreement between pure forward sampler vs.
 * forward+end point sampler.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class TestForwardAndEndPointSamplers
{
  public static void main(String [] args)
  {
    new TestForwardAndEndPointSamplers().test();
  }
  
  @Test
  public void test()
  {
    SimpleRateMatrix k80 = RateMatrices.kimura1980();
    Random rand = new Random(1);
    final double T = 3.0;
    
    PathStatistics pathStat1 = new PathStatistics(k80.nStates());
    PathStatistics pathStat2 = new PathStatistics(k80.nStates());
    
    CTMC process = k80.getProcess();
    ForwardSampler fwdSampler   = new ForwardSampler(process);
    EndPointSampler postSampler = new EndPointSampler(process);
    double nIters = 1000000;
    for (int i = 0; i < nIters; i++)
    {
      Path current = new Path();
      fwdSampler.sample(rand, T, pathStat1, current);
      Path p2 = new Path();
      postSampler.sample(rand, current.firstState(), current.lastState(), T, pathStat2, p2);
    }
    System.out.println("Cache size: " + postSampler.cacheSize());
    
    SimpleMatrix m1 = pathStat1.getCountsAsSimpleMatrix().scale(1.0/nIters);
    SimpleMatrix m2 = pathStat2.getCountsAsSimpleMatrix().scale(1.0/nIters);
    System.out.println(EJMLUtils.toString(m1));
    System.out.println();
    System.out.println(EJMLUtils.toString(m2));
    EJMLUtils.checkIsClose(m1, m2, 0.05);
  }
}
