package conifer.ctmc;

import tutorialj.Tutorial;

/**
 * A continuous time Markov chain. The main functionalities consists
 * in computing a marginal transition probability and a stationary distibution
 * (see below).
 * 
 * This implementation is based on caching the eigendecomposition
 * of the provided rate matrix.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class CTMC
{
  private final RateMatrix rateMatrix;
 
  /**
   * Note: if the RateMatrix is changed in place,
   * these changes will not be mirrored by this class.
   * 
   * It should be recreated each time a likelihood 
   * calculation is performed.
   * @param rateMatrix
   */
  public CTMC(RateMatrix rateMatrix)
  {
    this.rateMatrix = rateMatrix;
  }
  
  /**
   * 
   * Fill in ``marginalTransitionProbability()``.
   * 
   * Use the diagonalization method covered in class, using 
   * the eigen-decomposition functionalities provided by EJML.
   * 
   * See EJMLUtils.java in the bayonet project (v. 2.0.3).
   */
  @Tutorial(showSource = false, showLink = true)
  public double [][] marginalTransitionProbability(double branchLength)
  {
    throw new RuntimeException();
 
  }

  /**
   * Fill in ``stationaryDistribution()``.
   * 
   * Show how you can reduce the problem of finding the stationary distribution
   * of a CTMC to a eigenvalue problem. Again, solve this problem using
   * EJML.
   */
  @Tutorial(showSource = false, showLink = true)
  public double [] stationaryDistribution()
  {
    throw new RuntimeException();

  }
  

}
