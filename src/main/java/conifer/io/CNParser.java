package conifer.io;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import briefj.BriefIO;
import briefj.Indexer;
import conifer.models.CNPair;
import conifer.models.CNSpecies;

import com.google.common.collect.Maps;

import conifer.TreeNode;

/**
 * Create TreeNode to CharSequence map
 * 
 * @author Sohrab Salehi (sohrab.salehi@gmail.com)
 *
 */
public class CNParser {

	// Expecting the CSV file to have the following header:
	// species, siteName, ra, rA, cluster
	// and an optional header	
	public static Set<CNSpecies> readCNPairs(File file) 
	{
		Set<CNSpecies> cnSpecies = new LinkedHashSet<CNSpecies>();
		String WHITE_SPACE = "(\\s+$|^\\s+)";

		// CSVReader reader = new CSVReader(new FileReader(f));
		int lineNumber = 0;

		String currentTaxonName = null;
		List<CNPair> currentPairs = null;

		String clusterID = "";
		String speciesName = null;

		for (String line : BriefIO.readLines(file)) {
			if (line.length() > 0) {
				// nextLine[] is an array of values from the line
				// TODO: make it more general
				String[] elements = line.split(",");

				// skip header & comments & malformed elements
				if (line.charAt(0) == '#' || elements.length < 3 || isHeader(elements)) {
					System.out.println("This line was skipped: " + line);
					continue;
				} 
				
				// parse the line
				// TODO: match column names
				speciesName = elements[0].replaceAll(WHITE_SPACE, "");
				CNPair currentPair = CNPair.withCounts(
						Integer.parseInt(elements[2].replaceAll(WHITE_SPACE, "")),
						Integer.parseInt(elements[3].replaceAll(WHITE_SPACE, "")));

				// get cluster_id if it's provided
				if (elements.length > 4) {
					clusterID = elements[4];
				}

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

					// add to the new representation 
					cnSpecies.add(new CNSpecies(currentPairs, clusterID, speciesName));

					// add the new item
					currentTaxonName = speciesName;
					currentPairs = new ArrayList<CNPair>();
					currentPairs.add(currentPair);
				}
			}
			lineNumber++;
		}

		// put the last list too
		if (currentPairs != null) { 
			cnSpecies.add(new CNSpecies(currentPairs, clusterID, speciesName));
		}

		return cnSpecies;
	}

//	public static LinkedHashMap<TreeNode, CharSequence> readCN(File f) {
//		LinkedHashMap<TreeNode, CharSequence> result = Maps.newLinkedHashMap();
//		LinkedHashMap<TreeNode, List<CNPair>> map2 = CNParser.readCNPairs(f);
//
//		for (TreeNode t : map2.keySet()) {
//			result.put(t, getString(map2.get(t)));
//		}
//
//		return result;
//	}

	public static String getString(List<CNPair> list) {
		StringBuilder current = new StringBuilder();
		for (CNPair p : list) {
			current.append(p.toString());
		}
		return current.toString();
	}

	// TODO: parse header and use appropriate indexes
	public static boolean isHeader(String[] s) {
		if (isNumber(s[2]) && isNumber(s[3])) return false;
		return true;
	}
	
	// return true if it's a number
	public static boolean isNumber(String s) {	
		return  s.matches("-?\\d+(\\.\\d+)?");
	}

	public static LinkedHashMap<TreeNode, List<CNPair>> getNodeMap(Set<CNSpecies> cnSpecies) {
		LinkedHashMap<TreeNode, List<CNPair>> map = Maps.newLinkedHashMap();
		for (CNSpecies s:cnSpecies) {
			map.put(TreeNode.withLabel(s.getSpeciesName()), s.getCnPairs());
		}
		return map;
	}

	public static Map<String, Map<TreeNode, List<CNPair>>> getTreeNodeClusters(Set<CNSpecies> cnSpecies) {
		Map<String, Map<TreeNode, List<CNPair>>> result = new LinkedHashMap<String, Map<TreeNode, List<CNPair>>>();

		for (CNSpecies s:cnSpecies) {
			String clusterID = s.getClusterID();
			Map<TreeNode, List<CNPair>> temp = result.get(clusterID);
			if (temp == null) {
				temp = new LinkedHashMap<TreeNode, List<CNPair>>();
			} 

			temp.put(TreeNode.withLabel(s.getSpeciesName()), s.getCnPairs());
			result.put(clusterID, temp);
		}
		return result;
	}
	
	public static Map<String, List<CNSpecies>> getSpeciesClusters(Set<CNSpecies> cnSpecies) {
		Map<String, List<CNSpecies>> result = new LinkedHashMap<String, List<CNSpecies>>();

		for (CNSpecies s:cnSpecies) {
			String clusterID = s.getClusterID();
			List<CNSpecies> tempCluster = result.get(clusterID);
			if (tempCluster == null) {
				tempCluster = new ArrayList<CNSpecies>();
			} 

			tempCluster.add(s);
			result.put(clusterID, tempCluster);
		}
		return result;
	}


	public static void main(String[] args) throws IOException {
		// File("src/main/resources/conifer/sampleInput/testCopyNumber.txt");
		File file = new File("src/main/resources/conifer/sampleInput/testPatientData.txt");

		Set<CNSpecies> cnSpecies = CNParser.readCNPairs(file);

		// check out the other representation
		for (CNSpecies s:cnSpecies) {
			System.out.println(s.toString());
		}
		
		// check LinkedHashMap<TreeNode, List<CNPair>> representation
		LinkedHashMap<TreeNode, List<CNPair>> leaves = CNParser.getNodeMap(cnSpecies);
		for (Map.Entry<TreeNode, List<CNPair>> e : leaves.entrySet()) {
			for (int i = 0; i < e.getValue().size(); i++) {
				System.out.println(e.getKey().toString() + ": " + e.getValue().get(i).toString());
			}
		}
		
		// check cluster representation for CNSpecies
		Map<String, List<CNSpecies>> clusters = CNParser.getSpeciesClusters(cnSpecies);
		for (String clusterID:clusters.keySet()) {
			System.out.println("ClusterID:" + clusterID);
			List<CNSpecies> l = clusters.get(clusterID);
			System.out.println(Arrays.deepToString(l.toArray()));
		}
		
		// check cluster representation for TreeNode, List<CNPair>
		Map<String, Map<TreeNode, List<CNPair>>> clustersT = CNParser.getTreeNodeClusters(cnSpecies);
		for (String clusterID:clustersT.keySet()) {
			System.out.println("ClusterID:" + clusterID);
			Map<TreeNode, List<CNPair>> l = clustersT.get(clusterID);
			Collection<List<CNPair>> ll = l.values();
			@SuppressWarnings("unchecked")
			Object[] lll = ll.toArray();
			for (Object cnpairs: lll) {
				System.out.println(cnpairs.toString());
			}
		}

		
		/*
		// example of an indexer
		Indexer<CNPair> indexer = Indexers.CNPairIndexer();
		System.out.println(indexer.i2o(5));

		// rawString representation
		LinkedHashMap<TreeNode, CharSequence> l = CNParser.readCN(f);
		for (Map.Entry<TreeNode, CharSequence> e : l.entrySet()) {
			System.out.println(e.getKey().toString() + ": " + e.getValue());
		}
		
		 */
		Indexer<String> indexer = Indexers.copyNumberCTMCIndexer();
		System.out.println(indexer.i2o(5));
		System.out.println(indexer.objectsList().size());
	}
}
