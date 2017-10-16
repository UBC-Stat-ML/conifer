package conifer.ctmc;




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
  
  /**
   * A model for transitions from latent states to observations,
   * or null if the latent and observation spaces coincide.
   * @return
   */
  public RateMatrixToEmissionModel getEmissionModel();

  public double[][] getRateMatrix();

}
