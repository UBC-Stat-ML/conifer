package conifer.RandomUtils;

import java.util.Random;

public interface PriorOfWeights<P>{
	
	public double logDensity(double x);
	public double gradientOflogDensity(double x);
	public double generate(Random random, P parameterization);
}
