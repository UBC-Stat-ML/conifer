package conifer.io;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.internal.Lists;

import conifer.ctmc.expfam.CTMCState;
import conifer.ctmc.expfam.SerializedExpFamMixture;
import conifer.ctmc.expfam.SerializedExpFamMixture.BinaryFeature;
import conifer.ctmc.expfam.SerializedExpFamMixture.UnaryFeature;

public class FeatureFactory {
	
	public static boolean fullSupport=true;
	public static int nCategories=1;
	public static final Object partition = null;
	public static final List<String> DNA = Lists.newArrayList("A", "C", "G", "T");
	public static final List<String> PROTEIN = Lists.newArrayList("A", "R", "N", "D", "C", "Q", "E", "G", "H", "L", "I", "K", "M", "F", "P", "S", "T", "W", "Y", "V");
	public static final List<String> CODON = PhylogeneticObservationFactory.codonFactory().orderedSymbols;
	private static List<List<String>> supportEdges = null;
	public static Map<String, String> codonsToAminoAcid = AminoAcidAndCodonMap.AminoAcidFromAndToCodons().getRight();
	
	
	public static List<UnaryFeature> constructUnaryFeatures(List<String> stateSpace){
		
		List<UnaryFeature> unaryFeatures = Lists.newArrayList();
		int categoryIndex = 0;
		for(String element : stateSpace){
			String latent = element;
			CTMCState ctmcState = new CTMCState(categoryIndex, latent, partition);
			String stateKeys = "statio(".concat(element).concat(")");
			Map<String, Double> features = new HashMap<>();
			features.put(stateKeys, Double.valueOf(1.0));
			unaryFeatures.add(new UnaryFeature(ctmcState, features));		
		}
		
		return unaryFeatures;
		
	}
	
	// code related to DNA modelling
	public static SerializedExpFamMixture dnaGTR(){
		
		List<String> orderedLatents = DNA;
		List<UnaryFeature> unaryFeatures = constructUnaryFeatures(DNA);
		List<Map<String, String>> featureTemplate = constructBinaryFeatureTemplate("GTR", DNA);
		List<BinaryFeature> binaryFeatures = constructAminoAcidPairwiseBinaryFeatureCombiningFeatureTemplates(featureTemplate, DNA);
	    return new SerializedExpFamMixture(nCategories, orderedLatents, supportEdges, unaryFeatures, binaryFeatures, fullSupport);
	}

	
	public static SerializedExpFamMixture kimura1980(){
	
    List<String> orderedLatents = DNA;
    List<UnaryFeature> unaryFeatures = constructUnaryFeatures(DNA);
    List<BinaryFeature> binaryFeatures = Lists.newArrayList();
    CTMCState ctmcState0 = new CTMCState((int) 0, "A", partition);
    CTMCState ctmcState1 = new CTMCState((int) 0, "G", partition);
    CTMCState ctmcState2 = new CTMCState((int) 0, "C", partition);
    CTMCState ctmcState3 = new CTMCState((int) 0, "T", partition);
    Map<String, Double> transitionFeature = new HashMap<>();
    transitionFeature.put("isTransition", Double.valueOf(1.0));
   
    BinaryFeature binaryFeature0 = new BinaryFeature(ctmcState0, ctmcState1, transitionFeature);
    BinaryFeature binaryFeature1 = new BinaryFeature(ctmcState2, ctmcState3, transitionFeature);
    binaryFeatures.add(binaryFeature0);
    binaryFeatures.add(binaryFeature1);
 
    return new SerializedExpFamMixture(nCategories, orderedLatents, supportEdges, unaryFeatures, binaryFeatures, fullSupport);		
			
}
	// code related to Amino Acid model
	public static Map<String, String> mapAminoAcidToPolarity(){
		
		// Acidic stands for Acidic Polar, Basic stands for Basic Polar, Non stands for NonPolar, Polar stands for Polar
		List<String> polarity = Lists.newArrayList("Acidic", "Basic", "Non", "Polar");
		Map<String, List<String>> polarityToAminoAcid = new HashMap<>();
		polarityToAminoAcid.put("Acidic", Lists.newArrayList("D", "E"));
		polarityToAminoAcid.put("Basic", Lists.newArrayList("R", "H", "K"));
		polarityToAminoAcid.put("Non", Lists.newArrayList("A", "C", "G", "I", "L", "M", "F", "P", "W", "V"));
		polarityToAminoAcid.put("Polar", Lists.newArrayList("N", "Q", "S", "T", "Y"));
		Map<String, String> aminoAcidToPolarity = new HashMap<>();
		for(String ele: polarityToAminoAcid.keySet()){
			List<String> values = polarityToAminoAcid.get(ele);
			for(String value:values){
				aminoAcidToPolarity.put(value, ele);
			}
			
		}		
		return aminoAcidToPolarity;
	}
	
	
	public static Map<String, String> mapAminoAcidToSize(){
		
		Map<String, List<String>> sizeToAminoAcid = new HashMap<>();
		sizeToAminoAcid.put("Micro", Lists.newArrayList("A", "G", "S", "P", "C", "T", "N", "D"));
		sizeToAminoAcid.put("Big", Lists.newArrayList("R", "E", "Q", "H", "I", "L", "K", "M", "F", "W", "Y", "V"));
		Map<String, String> aminoAcidToSize = new HashMap<>();
		for(String ele:sizeToAminoAcid.keySet()){
			List<String> values = sizeToAminoAcid.get(ele);
			for(String value:values){
				aminoAcidToSize.put(value, ele);
			}
		}
	
		return aminoAcidToSize;
		
	}
	
	public static Map<String, String> mapStateSpaceCharacterToGTR(List<String> stateSpace){
		Map<String, String>  charToGTR = new LinkedHashMap<>();
		for(String ele: stateSpace){
			charToGTR.put(ele, ele);
		}
		return charToGTR;
	}
	
	
	public static List<BinaryFeature> constructAminoAcidPairwiseBinaryFeatureCombiningFeatureTemplates(List<Map<String, String>> featureTemplates, List<String> stateSpace){
		
		List<BinaryFeature> binaryFeatures = Lists.newArrayList();
		List<String> wholeState = Lists.newArrayList();
		wholeState.addAll(stateSpace);
		int categoryIndex = 0;
		
		while(wholeState.size()>0){
			String state0 = wholeState.get(0);
			CTMCState ctmcState0 = new CTMCState(categoryIndex, state0, partition);
			wholeState.remove(state0);
			
			// loop over all feature templates
			for(String state1:wholeState){
				CTMCState ctmcState1 = new CTMCState(categoryIndex, state1, partition);
				Map<String, Double> features = new LinkedHashMap<>();
				
				for(Map<String, String> featureTemplate: featureTemplates){
					String propertyOfState0 = featureTemplate.get(state0);
					String propertyOfState1 = featureTemplate.get(state1);
					String featureKeys = (propertyOfState0.compareTo(propertyOfState1)< 0)? propertyOfState0.concat(propertyOfState1): propertyOfState1.concat(propertyOfState0);
					features.put(featureKeys, Double.valueOf(1.0));
					}
				BinaryFeature binaryFeature = new BinaryFeature(ctmcState0, ctmcState1, features);
				binaryFeatures.add(binaryFeature);
				
			}
			
		}
		
		return binaryFeatures;	
	}
	
	
	public static List<Map<String, String>> constructBinaryFeatureTemplate(String featureTemplateNames, List<String> stateSpace){
	    List<Map<String, String>> result = Lists.newArrayList();
		if(featureTemplateNames.contentEquals("polarity")){
			result.add( mapAminoAcidToPolarity());
		}else if(featureTemplateNames.contentEquals("size")){
			result.add(mapAminoAcidToSize());
		}else if(featureTemplateNames.contentEquals("polarityPlusSize")){
			result.add( mapAminoAcidToPolarity());
			result.add(mapAminoAcidToSize());
		}else if(featureTemplateNames.contentEquals("polarityPlusSizePlusGTR")){
			
			result.add(mapAminoAcidToPolarity());
			result.add(mapAminoAcidToSize());
			result.add(mapStateSpaceCharacterToGTR(stateSpace));
		}else if(featureTemplateNames.contentEquals("GTR")){
			result.add(mapStateSpaceCharacterToGTR(stateSpace));
		}
			
		return result;
	}

	public static SerializedExpFamMixture proteinSimpleGTR(){
		List<String> orderedLatents = PROTEIN;
		List<UnaryFeature> unaryFeatures = constructUnaryFeatures(PROTEIN);
		List<Map<String, String>> featureTemplate = constructBinaryFeatureTemplate("GTR", PROTEIN);
		List<BinaryFeature> binaryFeatures = constructAminoAcidPairwiseBinaryFeatureCombiningFeatureTemplates(featureTemplate, PROTEIN);
		return new SerializedExpFamMixture(nCategories, orderedLatents, supportEdges, unaryFeatures, binaryFeatures, fullSupport);
	}


	public static SerializedExpFamMixture polarity(){
		List<String> orderedLatents = PROTEIN;
		List<UnaryFeature> unaryFeatures = constructUnaryFeatures(PROTEIN);
		List<Map<String, String>> featureTemplate = constructBinaryFeatureTemplate("polarity", PROTEIN);
		List<BinaryFeature> binaryFeatures = constructAminoAcidPairwiseBinaryFeatureCombiningFeatureTemplates(featureTemplate, PROTEIN);
		return new SerializedExpFamMixture(nCategories, orderedLatents, supportEdges, unaryFeatures, binaryFeatures, fullSupport);
		
	}
	
	public static SerializedExpFamMixture polaritySize(){
		List<String> orderedLatents = PROTEIN;
		List<UnaryFeature> unaryFeatures = constructUnaryFeatures(PROTEIN);
		List<Map<String,String>> featureTemplate = constructBinaryFeatureTemplate("polarityPlusSize", PROTEIN);
		List<BinaryFeature> binaryFeatures = constructAminoAcidPairwiseBinaryFeatureCombiningFeatureTemplates(featureTemplate, PROTEIN);
		return new SerializedExpFamMixture(nCategories, orderedLatents, supportEdges, unaryFeatures, binaryFeatures, fullSupport); 
	}
	
	public static SerializedExpFamMixture polaritySizeGTR(){
		List<String> orderedLatents = PROTEIN;
		List<UnaryFeature> unaryFeatures = constructUnaryFeatures(PROTEIN);
		List<Map<String, String>> featureTemplate = constructBinaryFeatureTemplate("polarityPlusSizePlusGTR", PROTEIN);
		List<BinaryFeature> binaryFeatures = constructAminoAcidPairwiseBinaryFeatureCombiningFeatureTemplates(featureTemplate, PROTEIN);
		return new SerializedExpFamMixture(nCategories, orderedLatents, supportEdges, unaryFeatures, binaryFeatures, fullSupport); 
	}
	
	// code related to construct features for codon evolution
	// features at the Amino Acid level
	
	
	
	
	
	
	
	
	
	
	
	// features at the DNA level
	
	
	
	
	
	
	
	 public static void main(String [] args)
	    {
	        SerializedExpFamMixture s = SerializedExpFamMixture.fromResource("/conifer/ctmc/expfam/dnaGTR-expfam.txt");
	        SerializedExpFamMixture s1 = dnaGTR();
	        // compare the content of s and s1 are the same or not
	        System.out.println(s.toString().contentEquals(s1.toString()));
	        
	        SerializedExpFamMixture polaritySize = SerializedExpFamMixture.fromResource("/conifer/ctmc/expfam/polaritySize-expfam.txt");
	        SerializedExpFamMixture polaritySize1 = polaritySize();
	        System.out.println(polaritySize.toString().contentEquals(polaritySize1.toString()));	        
	    }
	
	

}
