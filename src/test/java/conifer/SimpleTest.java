package conifer;


import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import blang.MCMCRunner;
import blang.annotations.DefineFactor;



/**
 * Testing the BlangSmallExample model.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class SimpleTest extends MCMCRunner
{
  
  
  @DefineFactor 
  Exponential<RateParameterization> prior = Exponential.newExponential();
  
  @DefineFactor
  Exponential<RateParameterization> l0 = Exponential.newExponential().withRate(prior.realization);
  
  @DefineFactor
  Exponential<RateParameterization> l1 = Exponential.newExponential().withRate(prior.realization);
  
  @DefineFactor
  Exponential<RateParameterization> l2 = Exponential.newExponential().withRate(prior.realization);
  
  @DefineFactor
  Exponential<RateParameterization> l3 = Exponential.newExponential().withRate(prior.realization);
  
  @DefineFactor
  Exponential<RateParameterization> l4 = Exponential.newExponential().withRate(prior.realization);
  
  @DefineFactor
  Exponential<RateParameterization> l5 = Exponential.newExponential().withRate(prior.realization);
  
  @DefineFactor
  Exponential<RateParameterization> l6 = Exponential.newExponential().withRate(prior.realization);
  
  @DefineFactor
  Exponential<RateParameterization> l7 = Exponential.newExponential().withRate(prior.realization);
  

  public static void main(String[] args)
  {
    SimpleTest test = new SimpleTest();
    test.factory.mcmcOptions.nMCMCSweeps = 10000000;
    test.run();
  }
  

}
