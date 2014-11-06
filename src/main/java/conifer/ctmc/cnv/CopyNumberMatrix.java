package conifer.ctmc.cnv;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.ejml.simple.SimpleMatrix;

import bayonet.math.EJMLUtils;
import blang.annotations.FactorArgument;
import blang.variables.RealVariable;
import static blang.variables.RealVariable.real;
import conifer.ctmc.CTMC;
import conifer.ctmc.CTMCParameters;
import conifer.ctmc.EigenCTMC;
import conifer.ctmc.RateMatrixToEmissionModel;
import conifer.ctmc.RateMatrixUtils;
import conifer.ctmc.SimpleRateMatrix;

/**
 * Corresponds to notation in basic_model.pdf
 * 
 * 
 * @author Sean Jewell (jewellsean@gmail.com)
 *
 */
public class CopyNumberMatrix implements CTMCParameters
{
	@FactorArgument
	public final RealVariable alpha = real(1); 

	@FactorArgument
	public final RealVariable beta = real(1); 

	@FactorArgument
	public final RealVariable gamma = real(1); 

	private final double DOLLO_EPSILON = 0.001;  
	
	
	// Break the overall Q representation into a series of 3 matrices
	private final SimpleMatrix increaseQ; 
	private final SimpleMatrix decreaseQ; 
	private final SimpleMatrix mutationQ;
	
	private final RateMatrixToEmissionModel emissionModel;

	public CopyNumberMatrix(double[][] increaseQ, 
			                double[][] decreaseQ,
			                double[][] mutationQ,
			                RateMatrixToEmissionModel emissionModel)
	{
		this.mutationQ = new SimpleMatrix(mutationQ); 
		this.increaseQ = new SimpleMatrix(increaseQ); 
		this.decreaseQ = new SimpleMatrix(decreaseQ); 
		this.emissionModel = emissionModel;
	}



	@Override
	public CTMC getProcess()
	{
		return new EigenCTMC(this.getRateMatrix());
	}

	@Override
	public RateMatrixToEmissionModel getEmissionModel()
	{
		return emissionModel;
	}

	/**
	 * Need to normalize the combined rate matrix 
	 * then return a normalized Q multiplied by DOLLO_EPSILON
	 * this is needed for Stochastic Dollo intrepretation cf. Alex's convergence result
	 */
	@Override
	public double[][] getRateMatrix()
	{
		double[][] rate = EJMLUtils.copyMatrixToArray(
				mutationQ.scale(gamma.getValue())
				.plus(increaseQ.scale(alpha.getValue()))
				.plus(decreaseQ.scale(beta.getValue())));
		RateMatrixUtils.normalize(rate); 
		return EJMLUtils.copyMatrixToArray((new SimpleMatrix(rate)).scale(DOLLO_EPSILON)); 
	}

	/* Syntactic Candy */ 
	/**
	 * @Warning emission matrix not currently implemented! 
	 * @param resourceURLs
	 * @return CopyNumberMatrix
	 */
	public static CopyNumberMatrix fromResources(Map<String, String> resourceURLs)
	{
		return new CopyNumberMatrix(
				getResourceFromLabel("mutationQ", resourceURLs),
				getResourceFromLabel("increaseQ", resourceURLs),
				getResourceFromLabel("decreaseQ", resourceURLs),
				null);
	}
	
	public static double[][] getResourceFromLabel(String label, 
			Map<String, String> resourceURLs)
	{
		return fromResource(resourceURLs.get(label));
	}

	public static double[][] fromResource(String resourceURL)
	{
		SimpleRateMatrix resourceMatrix = SimpleRateMatrix.fromResource(resourceURL);
		return resourceMatrix.getRateMatrix();
	}
	


	public static void main(String args [])
	{
		Set<String> labels = new HashSet<>();
		labels.add("mutationQ"); 
		labels.add("increaseQ");
		labels.add("decreaseQ");


		Map<String, String> resources = new HashMap<>();
		for(String lbl : labels)
			resources.put(lbl, "/conifer/ctmc/kimura1980.txt");

		CopyNumberMatrix cp = fromResources(resources);

		System.out.println((new SimpleMatrix(cp.getRateMatrix())).toString());

	}
}
