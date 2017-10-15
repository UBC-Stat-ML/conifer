package conifer.io;

import com.beust.jcommander.internal.Lists;
import java.util.*;

import org.apache.commons.lang3.tuple.Pair;

/**
 * AminoAcidToCompressedCodons() is a HashMap where the keys are the amino acids and the values are the ambiguous symbols of the codons
 * AminoAcidFromAndToCodons() is a bijection where the left of the Pair saves the amino acids as the keys and all its corresponding codons are saved for each key
 * the right of the Pair saves 61 key-values pairs where each key is a codon and the the value is the corresponding amino acids, the three stopping codons are not
 * included in the keys.
 *  @author Tingting Zhao (zhaott0416@gmail.com)
*/
public class AminoAcidAndCodonMap {
	
	public static final List<String> aminoAcids = Lists.newArrayList("A", "R", "N", "D", "C", "Q", "E", "G", "H", "L", "I", "K", "M", "F", "P", "S", "T", "W", "Y", "V");
	
	public static Map<String, Set<String>> AminoAcidToCompressedCodons(){
		
		Map<String, Set<String>> aminoAcidToCompressedCodes = new HashMap<>();
		aminoAcidToCompressedCodes.put("A", new HashSet<>(Lists.newArrayList("GCN")));
		aminoAcidToCompressedCodes.put("R", new HashSet<>(Lists.newArrayList("CGN", "MGR")));
		aminoAcidToCompressedCodes.put("N", new HashSet<>(Lists.newArrayList("AAY")));
		aminoAcidToCompressedCodes.put("D", new HashSet<>(Lists.newArrayList("GAY")));
		aminoAcidToCompressedCodes.put("C", new HashSet<>(Lists.newArrayList("TGY")));
		aminoAcidToCompressedCodes.put("Q", new HashSet<>(Lists.newArrayList("CAR")));
		aminoAcidToCompressedCodes.put("E", new HashSet<>(Lists.newArrayList("GAR")));
		aminoAcidToCompressedCodes.put("G", new HashSet<>(Lists.newArrayList("GGN")));
		aminoAcidToCompressedCodes.put("H", new HashSet<>(Lists.newArrayList("CAY")));
		aminoAcidToCompressedCodes.put("I", new HashSet<>(Lists.newArrayList("ATH")));
		aminoAcidToCompressedCodes.put("L", new HashSet<>(Lists.newArrayList("YTR", "CTN")));
		aminoAcidToCompressedCodes.put("K", new HashSet<>(Lists.newArrayList("AAR")));
		aminoAcidToCompressedCodes.put("M", new HashSet<>(Lists.newArrayList("ATG")));
		aminoAcidToCompressedCodes.put("F", new HashSet<>(Lists.newArrayList("TTY")));
		aminoAcidToCompressedCodes.put("P", new HashSet<>(Lists.newArrayList("CCN")));
		aminoAcidToCompressedCodes.put("S", new HashSet<>(Lists.newArrayList("TCN", "AGY")));
		aminoAcidToCompressedCodes.put("T", new HashSet<>(Lists.newArrayList("ACN")));
		aminoAcidToCompressedCodes.put("W", new HashSet<>(Lists.newArrayList("TGG")));
		aminoAcidToCompressedCodes.put("Y", new HashSet<>(Lists.newArrayList("TAY")));
		aminoAcidToCompressedCodes.put("V", new HashSet<>(Lists.newArrayList("GTN")));
		
		return aminoAcidToCompressedCodes;
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
				Set<String> uniqueCodonsFromCompressedCodon = new HashSet<>();
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
}
