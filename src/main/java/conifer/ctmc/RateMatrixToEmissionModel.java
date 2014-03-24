package conifer.ctmc;




public interface RateMatrixToEmissionModel 
{
  /**
   * In some cases, the states in the CTMC represented by getRateMatrix()
   * take value in a state that is different than the space of observation.
   * 
   * The matrixStatesToObservationProbabilities gives the mapping (probabilities)
   * from latent to observed.
   * 
   * Note that it is also often the case the the two spaces are the same. In this
   * case, this matrix can be used as a measurement error probability.
   * 
   * @return A matrix with rows being latent states and columns being observed.
   */
  public double [][] getMatrixStatesToObservationProbabilities();
  
//  public Indexer<String> getObservationsIndexer();
}
