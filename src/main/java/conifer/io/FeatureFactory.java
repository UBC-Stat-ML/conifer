package conifer.io;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.internal.Lists;

import blang.inits.DesignatedConstructor;
import blang.inits.Input;
import conifer.ctmc.expfam.CTMCState;
import conifer.ctmc.expfam.SerializedExpFamMixture;
import conifer.ctmc.expfam.SerializedExpFamMixture.BinaryFeature;
import conifer.ctmc.expfam.SerializedExpFamMixture.UnaryFeature;
import conifer.io.featurefactory.AminoAcidAndCodonMap;
import conifer.io.featurefactory.MapAminoAcidToPolarity;
import conifer.io.featurefactory.MapAminoAcidToSize;
import conifer.io.featurefactory.MapSingleStateToFeatures;
import conifer.io.featurefactory.MapStateSpaceCharacterToGTR;

public class FeatureFactory {
	
	public static boolean fullSupport=true;
	public static List<String> stateSpace = null;
	public static int nCategories=1;
	public static final Object partition = null;
	public static final List<String> DNA = PhylogeneticObservationFactory.nucleotidesFactory().orderedSymbols;
	public static final List<String> PROTEIN = PhylogeneticObservationFactory.proteinFactory().orderedSymbols;
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
		
		return SerializedExpFamMixture.fromResource("/conifer/ctmc/expfam/dnaGTR-expfam.txt");
	}

	
	public static SerializedExpFamMixture kimura1980(){
		return SerializedExpFamMixture.fromResource("/conifer/ctmc/expfam/kimura1980-expfam.txt");		
			
    }
	// code related to Amino Acid model
	public static Map<String, String> mapAminoAcidToPolarity() {
		
		MapSingleStateToFeatures featureClass = new MapAminoAcidToPolarity();
		return featureClass.getStatesAndFeatures();
	}
	
	
	public static Map<String, String> mapAminoAcidToSize(){
		MapSingleStateToFeatures featureClass = new MapAminoAcidToSize();
		return featureClass.getStatesAndFeatures();
		
	}
	
	@DesignatedConstructor
	public static void setStateSpace(@Input(formatDescription = "DNA, protein, codon or path to JSON spec") String description){
		String cleanedDescr = description.trim().toUpperCase();
	    if (cleanedDescr.equals("DNA"))
	      stateSpace = DNA;
	    else if (cleanedDescr.equals("PROTEIN"))
	      stateSpace = PROTEIN;
	    else if (cleanedDescr.equals("CODON"))
	      stateSpace = CODON;
	    }
	
	public static Map<String, String> mapStateSpaceCharacterToGTR(String description){
		
		MapSingleStateToFeatures featureClass = new MapStateSpaceCharacterToGTR();
		setStateSpace(description);
		featureClass.setStateSpace(stateSpace);
		return featureClass.getStatesAndFeatures();
	}
	
	
	public static List<BinaryFeature> constructPairwiseBinaryFeatureCombiningFeatureTemplates(List<Map<String, String>> featureTemplates, List<String> stateSpace){
		
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
	
	
	public static List<Map<String, String>> constructBinaryFeatureTemplate(String featureTemplateNames, String description){
		// description can be DNA, PROTEIN and CODON
	    List<Map<String, String>> result = Lists.newArrayList();
		if(featureTemplateNames.contentEquals("polarity")){
			result.add( mapAminoAcidToPolarity());
		}else if(featureTemplateNames.contentEquals("size")){
			result.add(mapAminoAcidToSize());
		}else if(featureTemplateNames.contentEquals("polaritySize")){
			result.add( mapAminoAcidToPolarity());
			result.add(mapAminoAcidToSize());
		}else if(featureTemplateNames.contentEquals("polaritySizeGTR")){
			
			result.add(mapAminoAcidToPolarity());
			result.add(mapAminoAcidToSize());
			setStateSpace(description);
			result.add(mapStateSpaceCharacterToGTR(description));
		}else if(featureTemplateNames.contentEquals("GTR")){
			setStateSpace(description);
			result.add(mapStateSpaceCharacterToGTR(description));
		}
			
		return result;
	}

	public static SerializedExpFamMixture proteinSimpleGTR(){
		List<String> orderedLatents = PROTEIN;
		List<UnaryFeature> unaryFeatures = constructUnaryFeatures(PROTEIN);
		List<Map<String, String>> featureTemplate = constructBinaryFeatureTemplate("GTR", "PROTEIN");
		List<BinaryFeature> binaryFeatures = constructPairwiseBinaryFeatureCombiningFeatureTemplates(featureTemplate, PROTEIN);
		return new SerializedExpFamMixture(nCategories, orderedLatents, supportEdges, unaryFeatures, binaryFeatures, fullSupport);
	}

	public static SerializedExpFamMixture constructEvolutionModel(String modelNames, List<String> stateSpace, String stateSpaceDescription){
		
		List<String> orderedLatents = stateSpace;
		List<UnaryFeature> unaryFeatures = constructUnaryFeatures(stateSpace);
		List<Map<String, String>> featureTemplate = constructBinaryFeatureTemplate(modelNames, stateSpaceDescription);
		List<BinaryFeature> binaryFeatures = constructPairwiseBinaryFeatureCombiningFeatureTemplates(featureTemplate, stateSpace);
		return new SerializedExpFamMixture(nCategories, orderedLatents, supportEdges, unaryFeatures, binaryFeatures, fullSupport);	
	}

	public static SerializedExpFamMixture polarity(){
		
		return constructEvolutionModel("polarity", PROTEIN, "PROTEIN");
		
	}
	
	public static SerializedExpFamMixture polaritySize(){
		
		return constructEvolutionModel("polaritySize", PROTEIN, "PROTEIN");
		
	}
	
	public static SerializedExpFamMixture polaritySizeGTR(){
		
		return constructEvolutionModel("polaritySizeGTR", PROTEIN, "PROTEIN");
		
	}
	
	
	// code related to construct features for codon evolution
	// construct codon level features, which is synonymous/non-synonymous feature
	public static Map<String, Map<String, String>> mapPairOfStatesToSynonymousFeature(String state0, String state1, List<String> stateSpace){
		// the default input for stateSpace should be codons
		String featureKeys = "isNonSynonymous";
		if(!stateSpace.contains(state0) || !stateSpace.contains(state1)){
			throw new RuntimeException( state0 + " or " + state1 + "does not belong to the statespace");
		}
		
		String aminoAcid0 = codonsToAminoAcid.get(state0);
		String aminoAcid1 = codonsToAminoAcid.get(state1);
		if(aminoAcid0.contentEquals(aminoAcid1)){
			featureKeys = "isSynonymous";
		}
		
		// return a data structure that can return state0, state1, and the corresponding featureKeys
		Map<String, Map<String, String>> result = new HashMap<>();
		Map<String, String> value = new HashMap<>();
		value.put(state1, featureKeys);
		result.put(state0, value);
		return result;
	} 
	
	public static Map<String, Map<String, String>> mapPairOfStatesToSynonymousFeature(String state0, String state1){
		
		return mapPairOfStatesToSynonymousFeature(state0, state1, CODON);
	}
	
	// features at the Amino Acid level
	
	public static Map<String, Map<String, String>> combineOneStateBijectionToPairOfStatesFeatureKeys(String state0, String state1, List<String> stateSpace, MapSingleStateToFeatures featureTemplateClass, String KeyOfState0, String KeyOfState1){
		
		featureTemplateClass.setStateSpace(stateSpace);
		Map<String, String> mapAllStates = featureTemplateClass.getStatesAndFeatures();
				
		String part1OfFeature = mapAllStates.get(state0);
		String part2OfFeature = mapAllStates.get(state1);
		
		// combine part1OfFeature and part2OfFeature according to their alphabetical order
		String featureKeys = (part1OfFeature.compareTo(part2OfFeature)< 0)? part1OfFeature.concat(part2OfFeature): part2OfFeature.concat(part1OfFeature);
		Map<String, Map<String, String>> result = new HashMap<>();
		Map<String, String> value = new HashMap<>();
		value.put(KeyOfState1, featureKeys);
		result.put(KeyOfState0, value);				
		return result;
	}
	
	public static Map<String, Map<String, String>> combineOneStateBijectionToPairOfStatesFeatureKeys(String state0, String state1, List<String> stateSpace, MapSingleStateToFeatures featureTemplateClass){
		
		return combineOneStateBijectionToPairOfStatesFeatureKeys(state0, state1, stateSpace, featureTemplateClass, state0, state1);
		
	}
	
	
	public static Map<String, Map<String, List<String>>> combineOneStateBijectionToPairOfStatesFeatureKeys(String state0, String state1, List<String> stateSpace, List<MapSingleStateToFeatures> featureTemplateClasses){
		
		Map<String, Map<String, List<String>>> result = new HashMap<>();
		List<String> values = Lists.newArrayList();
		for(MapSingleStateToFeatures element : featureTemplateClasses){
			Map<String, Map<String, String>> singleFeatureKey = (combineOneStateBijectionToPairOfStatesFeatureKeys(state0, state1, stateSpace, element));
			values.add(singleFeatureKey.get(state0).get(state1));
		}
		
		Map<String, List<String>> valuesOfKeyState0 = new HashMap<>();
		valuesOfKeyState0.put(state1, values);
		result.put(state0, valuesOfKeyState0);
		return result;
	}
	
	public static Map<String, Map<String, List<String>>> combineMultipleFeatureKeysOfTwoStates(String state0, String state1, Map<String, Map<String, List<String>>> listOfKeyFeatures, Map<String, Map<String, String>> newFeatureKeys){
		
		List<String> values = listOfKeyFeatures.get(state0).get(state1);
		values.add(newFeatureKeys.get(state0).get(state1));
		Map<String, List<String>> newValues = new HashMap<>();
		newValues.put(state1, values);
		listOfKeyFeatures.put(state0, newValues);	
		return listOfKeyFeatures;
	}
	
	// generate all amino acid level features for any two states
	public static Map<String, Map<String, String>> mapPairOfCodonsToPolarity(String state0, String state1, MapAminoAcidToPolarity mapAminoAcidToPolarity){
		// map the codon to the corresponding amino acid
		String aminoAcidOfState0 = codonsToAminoAcid.get(state0);
		String aminoAcidOfState1 = codonsToAminoAcid.get(state1);
		return combineOneStateBijectionToPairOfStatesFeatureKeys(aminoAcidOfState0, aminoAcidOfState1, PROTEIN, mapAminoAcidToPolarity, state0, state1);
	
	}
	
	public static Map<String, Map<String, String>> mapPairOfCodonsToSize(String state0, String state1, MapAminoAcidToSize mapAminoAcidToSize){
		// map the codon to the corresponding amino acid
		String aminoAcidOfState0 = codonsToAminoAcid.get(state0);
		String aminoAcidOfState1 = codonsToAminoAcid.get(state1);
		return combineOneStateBijectionToPairOfStatesFeatureKeys(aminoAcidOfState0, aminoAcidOfState1, PROTEIN, mapAminoAcidToSize, state0, state1);
	}
	
	public static Map<String, Map<String, String>> mapPairOfCodonsToGTR(String state0, String state1, MapStateSpaceCharacterToGTR mapStateSpaceCharacterToGTR){
		// map the codon to the corresponding amino acid
		// here the GTR feature, we refer to the GTR feature of the corresponding amino acid
		String aminoAcidOfState0 = codonsToAminoAcid.get(state0);
		String aminoAcidOfState1 = codonsToAminoAcid.get(state1);
		if(!aminoAcidOfState0.contentEquals(aminoAcidOfState1)){
			return combineOneStateBijectionToPairOfStatesFeatureKeys(aminoAcidOfState0, aminoAcidOfState1, PROTEIN, mapStateSpaceCharacterToGTR, state0, state1);			
		}else{
			Map<String, Map<String, String>> finalMap= new HashMap<>();
			Map<String,String> intemidiateMap = new HashMap<>();
			intemidiateMap.put(state1, "featureNonExistence");
			finalMap.put(state0, intemidiateMap);
			return finalMap;
		}
		
	
	}
		
	// generate all DNA level features
	public static List<String> splitCodonsIntoSingleCharacter(String codonState){
		List<String> result = Lists.newLinkedList();
		for(int i=0; i < codonState.length(); i++){
			Character element = codonState.charAt(i);
			String stringElement = String.valueOf(element);
			result.add(stringElement);
			
		}
		return result;
	}
	
	public static boolean checkIsTransition(String nucleotide0, String nucleotide1){
		boolean result = false;
		List<String> isTransitionPairs = Lists.newArrayList("AT", "TA", "CG", "GC");
		String twoStates = nucleotide0.concat(nucleotide1);
		if(isTransitionPairs.contains(twoStates)){
			result = true;
		}	
		return result;
	}
	
	
	public static Map<String, Map<String, List<String>>>  mapPairOfCodonsToTransitionPositionFeature(String state0, String state1){
		// When the state space is codon state, it consists of three nucleotide coming from A, C, G, T
		// We first split the state string into three letters
		
		if(state0.length()!=3 || state1.length()!=3){
			throw new RuntimeException("The state of codons does not consist of three nucleotide. The length is not equal to 3.");
		}
		
		List<String> threeLettersOfState0 = splitCodonsIntoSingleCharacter(state0);
		List<String> threeLettersOfState1 = splitCodonsIntoSingleCharacter(state1);
		// check each of the three positions and see if a transition happens
		List<String> featureStringNames = Lists.newLinkedList();
		featureStringNames.add(0, "isTransitionAt1stPosition");
		featureStringNames.add(1, "isTransitionAt2ndPosition");
		featureStringNames.add(2, "isTransitionAt3rdPosition");
		featureStringNames.add(3, "featureNonExistence");
		
		List<String> values = Lists.newArrayList();
		
		for(int i=0; i<state0.length(); i++){
			if(checkIsTransition(threeLettersOfState0.get(i), threeLettersOfState1.get(i))){
				// add feature names isTransition at the ith position
				values.add(featureStringNames.get(i));
			}		
		}
		if(values.size()==0){
			values.add(featureStringNames.get(3));
		}
		
		Map<String, Map<String, List<String>>> finalMap = new HashMap<>();
		Map<String, List<String>> intemediateMap = new HashMap<>();
		intemediateMap.put(state1, values);
		finalMap.put(state0, intemediateMap);
		return finalMap;
	}
	
	public static Map<String, Map<String, String>>  mapPairOfCodonsToSingleChangeFeature(String state0, String state1){
		
		// in this function, we assume only single change is allowed. 
		
		Map<String, Map<String, String>> result = new HashMap<>();
		Map<String, String> intemidiateMap = new HashMap<>();
		List<Boolean> labelPosition = Lists.newArrayList(false, false, false);
		
		// check the number of elements that are different
		// label the position when they are different
		int nTrues = 0;
		for(int i=0; i<state0.length(); i++){
			if(!String.valueOf(state0.charAt(i)).contentEquals(String.valueOf(state1.charAt(i)))){
				labelPosition.set(i, true);
				nTrues++;
			}			
		}
		
		List<String> featureStringNames = Lists.newArrayList("singleChangeAt1stPosition", "singleChangeAt2ndPosition", "singleChangeAt3rdPosition", "featureNonExistence");
		String value;
		
		// count the number of true element in a string
		if(Math.abs(nTrues-1) < 1e-6){
			// find the index of the element when it is true, if the index is i, then the (i+1)th position, there is a nucleotide change
			int j = labelPosition.indexOf(true);
			value = featureStringNames.get(j);
		
		}else {
			value = featureStringNames.get(3);
			
		}
		
		intemidiateMap.put(state1, value);
		result.put(state0, intemidiateMap);
		return result;
	}
	
	// create unary features and binary features using the feature keys for each pair of states to create codon models
	public static SerializedExpFamMixture codonUsingAminoAcidAndDNAFeatureModel(){
		List<String> orderedLatents = CODON;
		List<UnaryFeature> unaryFeatures = constructUnaryFeatures(CODON);
		List<BinaryFeature> binaryFeatures = Lists.newArrayList();
		
		List<String> wholeState = Lists.newArrayList();
		wholeState.addAll(CODON);
		int categoryIndex = 0;
		
		// combine all amino acid featureKeys into a list including polarity, size and GTR
		MapAminoAcidToPolarity mapAminoAcidToPolarity = new MapAminoAcidToPolarity();
		MapAminoAcidToSize mapAminoAcidToSize = new MapAminoAcidToSize();
		MapStateSpaceCharacterToGTR  mapStateSpaceCharacterToGTR = new MapStateSpaceCharacterToGTR();
		mapStateSpaceCharacterToGTR.stateSpace = CODON;
		
		// construct BinaryFeature for each pair of states
		while(wholeState.size()>0){
			
			String state0 = wholeState.get(0);
			CTMCState ctmcState0 = new CTMCState(categoryIndex, state0, partition);
			wholeState.remove(state0);
			
			for(String state1:wholeState){
				
				Map<String, Double> features = new LinkedHashMap<>();
				List<String> allFeatureMapKeys = Lists.newArrayList();
				CTMCState ctmcState1 = new CTMCState(categoryIndex, state1, partition);
				
				// put isSynonymous/isNonSynonymous feature in the template
				allFeatureMapKeys.add(mapPairOfStatesToSynonymousFeature(state0, state1, CODON).get(state0).get(state1));
				// put Amino Acid level feature in the template
				// add polarity, size and amino acid GTR features, combine all of them into a list
				allFeatureMapKeys.add(mapPairOfCodonsToPolarity(state0, state1, mapAminoAcidToPolarity).get(state0).get(state1));
				allFeatureMapKeys.add(mapPairOfCodonsToSize(state0, state1, mapAminoAcidToSize).get(state0).get(state1));
				allFeatureMapKeys.add(mapPairOfCodonsToGTR(state0, state1, mapStateSpaceCharacterToGTR).get(state0).get(state1));
				
			    // combine the existing features with features at DNA level
				allFeatureMapKeys.add(mapPairOfCodonsToSingleChangeFeature(state0, state1).get(state0).get(state1));
				allFeatureMapKeys.addAll(mapPairOfCodonsToTransitionPositionFeature(state0, state1).get(state0).get(state1));
				
				// Since if some feature is missing for the pair of state, we save the key string as "featureNonExistence",
				// as a result, we remove this string from the corresponding feature keys;
				while(allFeatureMapKeys.remove("featureNonExistence")){}
				
				
				// After getting all the features, we create the corresponding BinaryFeature according to the featureKeys
				for(String element: allFeatureMapKeys){
					features.put(element, Double.valueOf(1.0));
					}
				
				BinaryFeature binaryFeature = new BinaryFeature(ctmcState0, ctmcState1, features);
				binaryFeatures.add(binaryFeature);
			}
		}				
		
		return new SerializedExpFamMixture(nCategories, orderedLatents, supportEdges, unaryFeatures, binaryFeatures, fullSupport);
	}
		
	 public static void main(String [] args)
	    {
	
	        
	        SerializedExpFamMixture polaritySize = SerializedExpFamMixture.fromResource("/conifer/ctmc/expfam/polaritySize-expfam.txt");
	        SerializedExpFamMixture polaritySize1 = polaritySize();
	        System.out.println(polaritySize.toString());
	        System.out.println(polaritySize.toString().contentEquals(polaritySize1.toString()));	
	        
	        // check if the created codon model is correct in terms of the binary feature and unary feature
	        SerializedExpFamMixture codonModel = codonUsingAminoAcidAndDNAFeatureModel();
	        System.out.println(codonModel.toString());
	        
	    }
	
	

}
