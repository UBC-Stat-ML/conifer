package conifer.ctmc;

import briefj.Indexer;



/**
 * Note: this is an interface so that various parameterizations can be 
 * maintained.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public interface CTMCParameters
{
  public CTMC getProcess();
  
  // TODO: this should not be String: either a generic or a wrapper around a Map
  // note: sometimes you do just want indexer over strings
//  public Indexer<String> getRateMatrixIndexer();
  
  /**
   * A model for transitions from latent states to observations,
   * or null if the latent and observation spaces coincide.
   * @return
   */
  public RateMatrixToEmissionModel getEmissionModel();

  public double[][] getRateMatrix();

}
