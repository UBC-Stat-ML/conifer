package conifer.io.featurefactory;

import com.google.gson.Gson;

import briefj.BriefIO;
import conifer.io.PhylogeneticObservationFactory;
import com.google.common.collect.Sets;
import java.util.*;
import org.apache.commons.lang3.tuple.Pair;
import com.google.gson.reflect.TypeToken;

/**
 * AminoAcidToCompressedCodons() is a HashMap where the keys are the amino acids and the values are the ambiguous symbols of the codons
 * AminoAcidFromAndToCodons() is a bijection where the left of the Pair saves the amino acids as the keys and all its corresponding codons are saved for each key
 * the right of the Pair saves 61 key-values pairs where each key is a codon and the the value is the corresponding amino acids, the three stopping codons are not
 * included in the keys.
 *  @author Tingting Zhao (zhaott0416@gmail.com)
*/
public class AminoAcidAndCodonMap {
	
	public static final List<String> aminoAcids = PhylogeneticObservationFactory.proteinFactory().orderedSymbols;
	public static Map<String, Set<String>> aminoAcidToCompressedCodon = new HashMap<String, Set<String>>();
	 
	private static Map<String, Set<String>> _aminoAcidToCompressedCodons = null;
	
	private static Map<String, Set<String>> fromResource(String resourceURL)
	{
	    String jsonString = BriefIO.resourceToString(resourceURL, AminoAcidAndCodonMap.class); 
	    return fromJSONString(jsonString);
	    }
	  
	private static Map<String, Set<String>> fromJSONString(String jsonString)
	{
		return new Gson().fromJson(jsonString, new TypeToken<Map<String, Set<String>>>(){}.getType());
	    }
	
	public static Map<String, Set<String>> AminoAcidToCompressedCodons(){
		
		if (_aminoAcidToCompressedCodons == null)
			_aminoAcidToCompressedCodons = fromResource("/conifer/io/AminoAcidToCompressedCodons.txt");
		    return _aminoAcidToCompressedCodons;
		
	}
	
	public static Pair<Map<String, Set<String>>,Map<String, String>> AminoAcidFromAndToCodons(){
		Map<String, Set<String>> aminoAcidToCodons = new HashMap<>();
		Map<String, String> codonsToAminoAcids = new HashMap<>();
		
		Map<String, Set<String>> aminoAcidToCompressedCodon = AminoAcidToCompressedCodons();
		Map<String, Set<String>> ambiguousSymbols = PhylogeneticObservationFactory.codonFactory().ambiguousSymbols;
		Set<String> allCompressedCodonsKeySet = ambiguousSymbols.keySet();
		
		for(String ele:aminoAcids){
			
			// from compressedCodons, we get the exact codons, given a key from amino acid, we get all possible codons that can code this amino acid 
			Set<String> compressedCodons = aminoAcidToCompressedCodon.get(ele);
					
			// Among compressedCodons, ATG coded for amino acid M, TGG coded for W are not compressed code
			if(allCompressedCodonsKeySet.containsAll(compressedCodons)){
				// from key "compressedCodons" to get all original possible codons
				// the compressedCodons is a set, we need to loop over all of its elements
				Set<String> uniqueCodonsFromCompressedCodon = Sets.newHashSet();
			    for(String compressedCodon : compressedCodons){
					uniqueCodonsFromCompressedCodon.addAll( ambiguousSymbols.get(compressedCodon));
				}
			    // Key: amino acids, Values: uncompressed codons
			    aminoAcidToCodons.put(ele, uniqueCodonsFromCompressedCodon);
			    
			    // Key: uncompressed codons, Values: amino acids
				for(String codonEle : uniqueCodonsFromCompressedCodon){
					codonsToAminoAcids.put(codonEle, ele);
				}				
			}else{
				if(compressedCodons.size()>1){
					throw new RuntimeException("The number of codons which codes for M and W can't be bigger than one");
				}
				if(!compressedCodons.contains("ATG") && !compressedCodons.contains("TGG")){
					
					throw new RuntimeException("The codons coded for M and W can only be ATG and TGG");
				}
				
				aminoAcidToCodons.put(ele, compressedCodons);
				Iterator iter = compressedCodons.iterator();
				String first = (String) iter.next();
				
				if(compressedCodons.contains("ATG")||compressedCodons.contains("TGG"))
				codonsToAminoAcids.put(first, ele);			
			}			
			
		}
		// return a pair of result
		return Pair.of(aminoAcidToCodons, codonsToAminoAcids);
	}
	
	
	
	
	public static void main(String [] args){
		Map<String, Set<String>> aminoAcidToCompressedCodes = AminoAcidFromAndToCodons().getLeft();
	    System.out.println(aminoAcidToCompressedCodes.get("A").toString());
		String gson = new Gson().toJson(aminoAcidToCompressedCodes);
		System.out.println("Using reading json file method");
		System.out.println(JsonStringUtil.toPrettyFormat(gson));
		
	}
}


