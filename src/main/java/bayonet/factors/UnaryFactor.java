package bayonet.factors;



public interface UnaryFactor<V>
{
  public V connectedVariable();
  
  public double logNormalization();
}
