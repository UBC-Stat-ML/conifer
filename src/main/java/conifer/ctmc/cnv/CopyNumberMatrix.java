package conifer.ctmc.cnv;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.ejml.simple.SimpleMatrix;

import bayonet.math.EJMLUtils;
import blang.annotations.FactorArgument;
import blang.variables.ProbabilitySimplex;
import blang.variables.RealVariable;
import blang.variables.RealVector;
import briefj.opt.Option;
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

    // Ordered as: mutation, increase, decrease
    @FactorArgument
	public final RealVariable alpha = real(0.5);    


	@Option(gloss="Dollo precision parameter")
	public final double DOLLO_EPSILON = 10^-40;  
	
	
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
		double[] stationary = new double[increaseQ.numRows()];
		stationary[2] = 1; // this is the normal state
	    return new EigenCTMC(this.getRateMatrix(), stationary);
	}

	@Override
	public RateMatrixToEmissionModel getEmissionModel()
	{
		return emissionModel;
	}

	/**
	 * Need to normalize the combined rate matrix 
	 * then return a normalized Q multiplied by DOLLO_EPSILON
	 * this is needed for Stochastic Dollo interpretation cf. Alex's convergence result
	 */
	@Override
	public double[][] getRateMatrix()
	{
	    double param = alpha.getValue(); 
	    double[][] rate = EJMLUtils.copyMatrixToArray(
				mutationQ.plus(increaseQ.scale(param))
				.plus(decreaseQ.scale((1-param))).scale(DOLLO_EPSILON));
		RateMatrixUtils.fillRateMatrixDiagonalEntries(rate);
		return rate; 
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
	
	public static CopyNumberMatrix matrixOfSize(int size) 
	{	
		GenerateCNMatrices cn = new GenerateCNMatrices(size);
		cn.ensureInitalize();
		cn.saveMatricesAsJSON();
		return new CopyNumberMatrix(cn.getIncreaseQ(), cn.getDecreaseQ(), cn.getMutationQ(), null);
		/*
		Set<String> labels = new HashSet<>();
		labels.add("mutationQ"); 
		labels.add("increaseQ");
		labels.add("decreaseQ");

		Map<String, String> resources = new HashMap<>();
		for(String lbl : labels)
			resources.put(lbl, GenerateCNMatrices.getPathForMatrix(lbl, size));

		return fromResources(resources);
		*/
	}

	public static void main(String args [])
	{
		int size = 30;
		/*
		
		Set<String> labels = new HashSet<>();
		labels.add("mutationQ"); 
		labels.add("increaseQ");
		labels.add("decreaseQ");

		Map<String, String> resources = new HashMap<>();
		for(String lbl : labels) {
			resources.put(lbl, GenerateCNMatrices.getPathForMatrix(lbl, size));
		}
		CopyNumberMatrix cp = fromResources(resources);
		double[][] rateMatrix = cp.getRateMatrix();
		System.out.println(Arrays.deepToString(rateMatrix));
		System.out.println((new SimpleMatrix(rateMatrix)).toString());
		*/
		
		CopyNumberMatrix cp = CopyNumberMatrix.matrixOfSize(size);
		double[][] rateMatrix = cp.getRateMatrix();
		
		// check if it's a valid rate matrix:
		try {
			RateMatrixUtils.checkValidRateMatrix(rateMatrix);
			cp.getProcess();
		} catch(Exception e) {
			System.err.println(e);
			System.err.println(e.getStackTrace());
		}
	}
	
}
