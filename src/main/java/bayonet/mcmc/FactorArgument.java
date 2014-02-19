package bayonet.mcmc;


/**
 * Indicates a variable connected to a factor.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public @interface FactorArgument 
{
  /**
   * Whether this argument is part of the variable being
   * defined.
   * 
   * In other words, if X | Y ~ f(Y),
   * then both X and Y will be arguments of f, but only
   * X is being defined by f
   * 
   * @return
   */
  public boolean defines() default false;
}
