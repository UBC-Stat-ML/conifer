package conifer.io;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import briefj.BriefIO;
import briefj.Indexer;
import conifer.models.CNPair;

import com.google.common.collect.Maps;

import conifer.TreeNode;


public class CNParser {

	public static LinkedHashMap<TreeNode, List<CNPair>> readCN(File f) throws NumberFormatException, IOException {

/**
 * Create TreeNode to CharSequence map
 * @author Sohrab Salehi (sohrab.salehi@gmail.com)
 *
 */
public class CNParser {

	public static LinkedHashMap<TreeNode, List<CNPair>> readCNPairs(File f) {

		LinkedHashMap<TreeNode, List<CNPair>> map = Maps.newLinkedHashMap();
		
		/**
		 * We allow a line of comment, a line for header 
		 */
		
		String WHITE_SPACE = "(\\s+$|^\\s+)";
		
		//CSVReader reader = new CSVReader(new FileReader(f));
		int lineNumber = 0;
		
		String currentTaxonName = null;
		List<CNPair> currentPairs = null;
		
		for (String line : BriefIO.readLines(f))
		{
			if (line.length() > 0) {
				// nextLine[] is an array of values from the line
				
				// TODO: make it more general
				String[] elements = line.split(",");
				
				if (line.charAt(0) == '#' || lineNumber < 2) {
					// header & comment
					//System.out.println(Arrays.toString(line));
				} else {
					// TODO: add the site
					
					if (elements.length > 3) {
						
						// parse the line
						String speciesName = elements[0].replaceAll(WHITE_SPACE, "");
						CNPair currentPair = CNPair.withCounts(Integer.parseInt(elements[2].replaceAll(WHITE_SPACE, "")), 
								Integer.parseInt(elements[3].replaceAll(WHITE_SPACE, "")));
						
						// the first line
						if (currentTaxonName == null) {
							currentTaxonName = speciesName;
							currentPairs = new ArrayList<CNPair>();
							currentPairs.add(currentPair);
						} else if (currentTaxonName.equals(speciesName)) {
							// still the same species, continue adding the pairs
							currentPairs.add(currentPair);
						} else {
							// new species, put the last one, and renew
							TreeNode currentTaxon = TreeNode.withLabel(currentTaxonName);
							map.put(currentTaxon, currentPairs);
							
							// add the new item
							currentTaxonName = speciesName;
							currentPairs = new ArrayList<CNPair>();
							currentPairs.add(currentPair);
						}
					}
				}
			}
			lineNumber++;
		}
		
		// put the last list too
		if (currentPairs != null) {
			map.put(TreeNode.withLabel(currentTaxonName), currentPairs);
		}
		
		return map;
	}
	
	
	
	public static LinkedHashMap<TreeNode, CharSequence> readCN(File f)  {
		LinkedHashMap<TreeNode, CharSequence> result = Maps.newLinkedHashMap();
		LinkedHashMap<TreeNode, List<CNPair>> map2 = CNParser.readCNPairs(f);
		
		for (TreeNode t : map2.keySet()) {
			result.put(t, getString(map2.get(t)));
		}
		
		return result;
	}
	
	public static String getString(List<CNPair> list) {
		StringBuilder current = new StringBuilder();
		for (CNPair p : list) {
			current.append(p.toRawString());
		}
		return current.toString();
	}
	
	public static void main(String[] args) throws IOException {
		//File f = new File("src/main/resources/conifer/sampleInput/testCopyNumber.txt");
		File f = new File("src/main/resources/conifer/sampleInput/testPatientData.txt");
		LinkedHashMap<TreeNode, List<CNPair>> leaves = CNParser.readCNPairs(f);
		
		
		for (Map.Entry<TreeNode, List<CNPair>> e: leaves.entrySet()) {
			for (int i = 0; i < e.getValue().size(); i++) {
				System.out.println(e.getKey().toString() + ": " +  e.getValue().get(i).toString());
			}
		}
		
		// example of an indexer
		Indexer<CNPair> indexer = Indexers.CNPairIndexer();
		System.out.println(indexer.i2o(5));

		
		// rawString representation
		LinkedHashMap<TreeNode, CharSequence> l = CNParser.readCN(f);
		for (Map.Entry<TreeNode, CharSequence> e: l.entrySet()) {
			System.out.println(e.getKey().toString() + ": " +  e.getValue());
		}
		
	}
}
