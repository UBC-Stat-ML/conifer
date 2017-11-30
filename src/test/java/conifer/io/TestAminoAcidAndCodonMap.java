package conifer.io;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Assert;
import org.junit.Test;

import com.beust.jcommander.internal.Lists;

import conifer.io.featurefactory.AminoAcidAndCodonMap;


public class TestAminoAcidAndCodonMap
{
 
  @Test
  public void testAminoAcidToCompressedCodon()
  {
	  
	  Map<String, Set<String>> compressedCodons = AminoAcidAndCodonMap.AminoAcidToCompressedCodons();
	  // randomly choose some amino acids and see if the key and values are correct
	  Set<String> rValues = compressedCodons.get("R");
	  Set<String> lValues = compressedCodons.get("L");
	  Set<String> sValues = compressedCodons.get("S");
	  
	  Assert.assertEquals(rValues.size(), 2);
	  Assert.assertEquals(lValues.size(), 2);
	  Assert.assertEquals(sValues.size(), 2);
	  
	  Assert.assertTrue(rValues.contains("CGN"));
	  Assert.assertTrue(rValues.contains("MGR"));
	  Assert.assertFalse(rValues.contains("AAY"));
	  
	  Assert.assertTrue(lValues.contains("YTR"));
	  Assert.assertTrue(lValues.contains("CTN"));
	  
	  Assert.assertTrue(sValues.contains("TCN"));
	  Assert.assertTrue(sValues.contains("AGY"));
	  
	  Assert.assertTrue(compressedCodons.get("I").contains("ATH"));
	  
  }
  
  @Test
  public void testAminoAcidFromAndToCodons()
  {
	  Pair<Map<String, Set<String>>,Map<String, String>> result = AminoAcidAndCodonMap.AminoAcidFromAndToCodons();
	  Map<String, Set<String>> aminoAcidToAllCodons = result.getLeft();
	  Map<String, String> codonsToAminoAcids = result.getRight();
	  
	  // check if the codons for all amino acids in terms of the key value pairs is correct
	  // test amino acid as keys and codons as values
	  Assert.assertEquals(aminoAcidToAllCodons.get("L").size(), 6);
	  Assert.assertTrue(aminoAcidToAllCodons.get("L").containsAll(new HashSet<>(Lists.newArrayList("TTA", "TTG", "CTA", "CTG", "CTC", "CTT"))));
	  
	  Assert.assertEquals(aminoAcidToAllCodons.get("I").size(), 3);
	  Assert.assertTrue(aminoAcidToAllCodons.get("I").containsAll(new HashSet<>(Lists.newArrayList("ATT", "ATC", "ATA"))));
	  
	  Assert.assertEquals(aminoAcidToAllCodons.get("M").size(), 1);
	  Assert.assertTrue(aminoAcidToAllCodons.get("M").containsAll(new HashSet<>(Lists.newArrayList("ATG"))));
	  
	  Assert.assertEquals(aminoAcidToAllCodons.get("V").size(), 4);
	  Assert.assertTrue(aminoAcidToAllCodons.get("V").containsAll(new HashSet<>(Lists.newArrayList("GTA", "GTC", "GTG", "GTT"))));
	  
	  Assert.assertEquals(aminoAcidToAllCodons.get("S").size(), 6);
	  Assert.assertTrue(aminoAcidToAllCodons.get("S").containsAll(new HashSet<>(Lists.newArrayList("TCA", "TCC", "TCG", "TCT", "AGC", "AGT"))));
	  
      Assert.assertEquals(aminoAcidToAllCodons.get("P").size(), 4);
      Assert.assertTrue(aminoAcidToAllCodons.get("P").containsAll(new HashSet<>(Lists.newArrayList("CCA", "CCC", "CCG", "CCT"))));
	  
	  Assert.assertEquals(aminoAcidToAllCodons.get("T").size(), 4);
	  Assert.assertTrue(aminoAcidToAllCodons.get("T").containsAll(new HashSet<>(Lists.newArrayList("ACA", "ACC", "ACG", "ACT"))));
	  
	  
	  Assert.assertEquals(aminoAcidToAllCodons.get("A").size(), 4);
	  Assert.assertTrue(aminoAcidToAllCodons.get("A").containsAll(new HashSet<>(Lists.newArrayList("GCA", "GCC", "GCG", "GCT"))));
	  
	  Assert.assertEquals(aminoAcidToAllCodons.get("G").size(), 4);
	  Assert.assertTrue(aminoAcidToAllCodons.get("G").containsAll(new HashSet<>(Lists.newArrayList("GGA", "GGC", "GGG", "GGT"))));
	  
	  Assert.assertEquals(aminoAcidToAllCodons.get("R").size(), 6);
	  Assert.assertTrue(aminoAcidToAllCodons.get("R").containsAll(new HashSet<>(Lists.newArrayList("CGA", "CGC", "CGG", "CGT", "AGA", "AGG"))));
	  
	  List<String> twoCodonsAminoAcidGroup = Lists.newArrayList("F", "Y", "H", "Q", "N", "K", "D", "E", "C");
	  for(String ele: twoCodonsAminoAcidGroup){
		  Assert.assertEquals(aminoAcidToAllCodons.get(ele).size(), 2);
	  }
	  
	  Assert.assertTrue(aminoAcidToAllCodons.get("F").containsAll(new HashSet<>(Lists.newArrayList("TTC", "TTT"))));
	  Assert.assertTrue(aminoAcidToAllCodons.get("Y").containsAll(new HashSet<>(Lists.newArrayList("TAT", "TAC"))));
	  Assert.assertTrue(aminoAcidToAllCodons.get("H").containsAll(new HashSet<>(Lists.newArrayList("CAC", "CAT"))));
	  Assert.assertTrue(aminoAcidToAllCodons.get("Q").containsAll(new HashSet<>(Lists.newArrayList("CAA", "CAG"))));
	  Assert.assertTrue(aminoAcidToAllCodons.get("N").containsAll(new HashSet<>(Lists.newArrayList("AAT", "AAC"))));
	  Assert.assertTrue(aminoAcidToAllCodons.get("K").containsAll(new HashSet<>(Lists.newArrayList("AAA", "AAG"))));
	  Assert.assertTrue(aminoAcidToAllCodons.get("D").containsAll(new HashSet<>(Lists.newArrayList("GAC", "GAT"))));
	  Assert.assertTrue(aminoAcidToAllCodons.get("E").containsAll(new HashSet<>(Lists.newArrayList("GAA", "GAG"))));
	  Assert.assertTrue(aminoAcidToAllCodons.get("C").containsAll(new HashSet<>(Lists.newArrayList("TGT", "TGC"))));
	  
	  Assert.assertEquals(aminoAcidToAllCodons.get("W").size(), 1);
	  Assert.assertTrue(aminoAcidToAllCodons.get("W").containsAll(new HashSet<>(Lists.newArrayList("TGG"))));
	  
	  // the values should not include the three stopping codons
	  Assert.assertFalse(aminoAcidToAllCodons.values().contains("TAA")|| aminoAcidToAllCodons.values().contains("TAA")||aminoAcidToAllCodons.values().contains("TGA"));
	  
	  // test codons as keys and amino acid as values
	  // Find the complement between the two sets
	  // The correctness of the code has been checked manually by printing out all key-value pairs by the following code
	  // System.out.println(Arrays.deepToString(codonsToAminoAcids.entrySet().toArray()));
	 
	  Assert.assertEquals(61, codonsToAminoAcids.entrySet().size());
	  Assert.assertTrue(codonsToAminoAcids.get("TGG").contentEquals("W"));
	  Assert.assertTrue(codonsToAminoAcids.get("ATG").contentEquals("M"));  
  }
  
  
}
