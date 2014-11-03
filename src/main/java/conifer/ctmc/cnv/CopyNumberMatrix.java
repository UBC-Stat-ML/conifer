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
 * 
 * @author Sean Jewell (jewellsean@gmail.com)
 *
 */
public class CopyNumberMatrix implements CTMCParameters
{
	@FactorArgument
	public final RealVariable increaseCNV = real(1); 

	@FactorArgument
	public final RealVariable decreaseCNV = real(1); 

	private final double mutationEpsilon = 0.001; 

	// Break the overall Q representation into a series of 3 matrices
	private final SimpleMatrix mutationQ;
	private final SimpleMatrix increaseQ; 
	private final SimpleMatrix decreaseQ; 

	private final RateMatrixToEmissionModel emissionModel;

	public CopyNumberMatrix(double[][] mutationQ, double[][] increaseQ, 
			double[][] decreaseQ,RateMatrixToEmissionModel emissionModel)
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
	 * then return a normalized Q
	 */
	@Override
	public double[][] getRateMatrix()
	{
		double[][] rate = EJMLUtils.copyMatrixToArray(
				mutationQ.scale(mutationEpsilon)
				.plus(increaseQ.scale(increaseCNV.getValue()))
				.plus(decreaseQ.scale(decreaseCNV.getValue())));
		RateMatrixUtils.normalize(rate); 
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
		//TODO: create the corresponding emission model.
		// Is s (BetaBinomial hyperparameter, always fixed?)
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
	
	/**
	 * @ Test method. Return some default values.
	 * @param cnEmissionModel
	 * @return
	 */
	public static CopyNumberMatrix testSingleCellCopyNumber(CopyNumberEmissionModel cnEmissionModel) {
		// TODO Auto-generated method stub
		Set<String> labels = new HashSet<>();
		labels.add("mutationQ"); 
		labels.add("increaseQ");
		labels.add("decreaseQ");

		Map<String, String> resources = new HashMap<>();
		for(String lbl : labels)
			resources.put(lbl, "/conifer/ctmc/kimura1980.txt");

		return new CopyNumberMatrix(
				getResourceFromLabel("mutationQ", resources),
				getResourceFromLabel("increaseQ", resources),
				getResourceFromLabel("decreaseQ", resources),
				cnEmissionModel);
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
