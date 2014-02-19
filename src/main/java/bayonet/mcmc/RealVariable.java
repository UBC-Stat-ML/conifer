package bayonet.mcmc;


/**
 * A simple real variable that can be modified in place. Used by
 * potentials to compute the log density of proposed/current states.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public interface RealVariable 
{
  /**
   * 
   * @param newValue
   */
  public void setValue(double newValue);
  
  /**
   * 
   * @return
   */
  public double getValue();
}