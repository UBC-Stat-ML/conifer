package conifer.ctmc.expfam;

import java.util.Random;

public interface WeightPrior<P>{
	
	public double logDensity(double x);
	public double gradientOflogDensity(double x);
	public double generate(Random random, P parameterization);
}
