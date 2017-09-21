package conifer.factors;

import java.util.Random;

public interface RealValuedDensity 
{
  public double logDensity(double point);
  public double sample(Random random);
}
