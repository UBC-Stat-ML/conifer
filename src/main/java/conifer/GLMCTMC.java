package conifer;
import static conifer.io.Indexers.dnaIndexer;
import static conifer.io.PhylogeneticObservationFactory.nucleotidesFactory;
import static conifer.io.PhylogeneticObservationFactory.proteinFactory;
import static conifer.io.PhylogeneticObservationFactory.codonFactory;
import static conifer.io.Indexers.proteinIndexer;

import java.util.Collections;
import java.util.Random;

import blang.core.RealConstant;
import blang.mcmc.SampledVariable;
import briefj.Indexer;
import conifer.EvoGLM;
import conifer.ctmc.RateMatrices;
import conifer.ctmc.RateMatrixToEmissionModel;
import conifer.ctmc.SimpleRateMatrix;
import conifer.ctmc.expfam.CTMCExpFam;
import conifer.ctmc.expfam.CTMCState;
import conifer.ctmc.expfam.CTMCStateSpace;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.ExpFamParameters;
import conifer.ctmc.expfam.SerializedExpFamMixture;
import conifer.io.PhylogeneticObservationFactory;
import conifer.models.DiscreteGammaMixture;
import conifer.models.MultiCategorySubstitutionModel;

public class GLMCTMC {
	
	// @SampledVariable(skipFactorsFromSampledModel = true)
	// public EvoGLM fullModel;
	
	public String modelNames;
	
	public String stateSpace;
	
	public GLMCTMC( String modelNames, String stateSpace) {
		// this.fullModel = fullModel;
		this.modelNames = modelNames;
		this.stateSpace = stateSpace;		
	}
	
	
	public SerializedExpFamMixture getSerialized(String modelNames, String stateSpace)
    {
		if (modelNames.equals("GTR") && stateSpace.trim().toUpperCase().equals("DNA"))
			return SerializedExpFamMixture.dnaGTR();
		else if (modelNames.equals("GTR") && stateSpace.trim().toUpperCase().equals("PROTEIN"))
			return SerializedExpFamMixture.proteinSimpleGTR();
		else if (modelNames.equals("polarity") && stateSpace.trim().toUpperCase().equals("PROTEIN"))
			return SerializedExpFamMixture.polarity();
		else if (modelNames.equals("polaritySize") && stateSpace.trim().toUpperCase().equals("PROTEIN"))
			return SerializedExpFamMixture.polaritySize();
		else if (modelNames.equals("polaritySizeGTR") && stateSpace.trim().toUpperCase().equals("PROTEIN"))
			return SerializedExpFamMixture.polaritySizeGTR();
		else if (modelNames.equals("kimura")&& stateSpace.trim().toUpperCase().equals("DNA"))
			return SerializedExpFamMixture.kimura1980();
		else {
			throw new RuntimeException("The name of the model is currently not supported");
		}
		//ToDo: add the codon case
		
		
    }

    public Indexer<String> getIndexer(String stateSpace)
    {
    		//ToDo: add the codon case
    		String cleanedDescr = stateSpace.trim().toUpperCase();
	    if (cleanedDescr.equals("DNA"))
	    		return dnaIndexer();
	    else if (cleanedDescr.equals("PROTEIN"))
	    		return proteinIndexer();
	    //else if (cleanedDescr.equals("CODON"))
	    		//return codonIndexer();
	    else {
	    		throw new RuntimeException("The name of the state space string is not supported");	    	
	    }	    
    }

    public PhylogeneticObservationFactory getFactory(String stateSpace)
    {
    		String cleanedDescr = stateSpace.trim().toUpperCase();
	    if (cleanedDescr.equals("DNA"))
	      return nucleotidesFactory();
	    else if (cleanedDescr.equals("PROTEIN"))
	      return proteinFactory();
	    else if (cleanedDescr.equals("CODON"))
	    	  return codonFactory();
	    else {
	    	  throw new RuntimeException("The name of the state space string is not supported");
	    }
    }
    
    
    
    
    public ExpFamMixture getMixture(String modelNames, String stateSpace) 
    {   	
    		SerializedExpFamMixture serialized = getSerialized(modelNames, stateSpace);
    		Indexer<String> observationIndexer = getIndexer(stateSpace);
    		boolean isNormalized = true;
    		ExpFamMixture mixture = ExpFamMixture.fromSerialized(serialized, observationIndexer, isNormalized);
    		return mixture;
    		    	
    }
    

    public SimpleRateMatrix getRateMtx(String modelNames, String stateSpace)
    {
    		// not sure of this, ask Alex
    		ExpFamMixture mixture = getMixture(modelNames, stateSpace);
    		ExpFamParameters param = mixture.parameters;
    		//param.setVector(values);
    		double [][] rateMtx = param.getRateMatrix(0);
    		SimpleRateMatrix result = new SimpleRateMatrix(rateMtx,mixture.emissionModel);
    		return result;
    	
    }
    
    @SuppressWarnings({ "rawtypes", "unchecked" })
    public MultiCategorySubstitutionModel<ExpFamMixture> getSubstitutionModel(String modelNames, String stateSpace,  int nSites)
    {
    		
    		ExpFamMixture mixture = getMixture(modelNames, stateSpace);
    		MultiCategorySubstitutionModel<ExpFamMixture> model = new MultiCategorySubstitutionModel(mixture, nSites);    
    		return model;   	
    }	

}
