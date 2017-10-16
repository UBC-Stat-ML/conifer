package conifer.ctmc;



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
public interface CTMC
{
  public double [][] marginalTransitionProbability(double branchLength);
  public double [] stationaryDistribution();
  public double [][] getRateMatrix();
}
