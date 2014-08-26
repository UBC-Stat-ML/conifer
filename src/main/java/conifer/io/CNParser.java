package conifer.io;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import briefj.BriefIO;
import conifer.models.CNPair;

import com.google.common.collect.Maps;

import conifer.TreeNode;

public class CNParser {

	public static LinkedHashMap<TreeNode, List<CNPair>> readCN(File f) throws NumberFormatException, IOException {
		LinkedHashMap<TreeNode, List<CNPair>> map = Maps.newLinkedHashMap();
		
		/**
		 * We allow a line of comment, a line for header 
		 */
		
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
						String speciesName = elements[0].replaceAll("(\\s+$|^\\s+)", "");
						CNPair currentPair = CNPair.withCounts(Integer.parseInt(elements[2].replaceAll("(\\s+$|^\\s+)", "")), 
								Integer.parseInt(elements[3].replaceAll("(\\s+$|^\\s+)", "")));
						
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
	
	
	public static void main(String[] args) throws IOException {
		File f = new File("src/main/resources/conifer/sampleInput/testCopyNumber.txt");
		LinkedHashMap<TreeNode, List<CNPair>> leaves = CNParser.readCN(f);
		
		for (Map.Entry<TreeNode, List<CNPair>> e: leaves.entrySet()) {
			for (int i = 0; i < e.getValue().size(); i++) {
				System.out.println(e.getKey().toString() + ": " +  e.getValue().get(i).toString());
			}
		}
		
		
	}
}
