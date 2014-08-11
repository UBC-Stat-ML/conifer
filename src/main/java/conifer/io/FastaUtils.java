package conifer.io;

import java.io.File;
import java.util.Arrays;
import java.util.LinkedHashMap;

import briefj.BriefIO;

import com.google.common.collect.Maps;

import conifer.TreeNode;


/**
 * Reads sequences in the FASTA format.
 * 
 * The first line of blocks of lines starting with '>' will become the sequence label.
 * 
 * The contents of the lines in between will be concatenated (ignoring blank sequence) and will
 * become the sequence associated with the the last block of comments.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class FastaUtils
{
	public static LinkedHashMap<TreeNode, CharSequence> readFasta(File f)
	{
		LinkedHashMap<TreeNode,CharSequence> map = Maps.newLinkedHashMap();
		StringBuilder current = null;
		for (String line : BriefIO.readLines(f))
		{
			if (line.length() > 0)
			{
				if (line.charAt(0) == '>' && 
						(current == null || current.length() > 0)) // this second condition will ignore all but the first line of comment blocks (usual convention)
				{
					current = new StringBuilder();
					TreeNode currentTaxon = TreeNode.withLabel(line.substring(1).replaceAll("(\\s+$|^\\s+)", ""));
					map.put(currentTaxon, current);
				}
				else
				{
					line = line.replaceAll("\\s+", "");
					if (current == null)
						throw new RuntimeException("Error in the currenlty picky implementation of the fasta reader: the first line did not start with '>'");
					current.append(line);
				}
			}
		}


		return map;
	}

	/**
	 * Currently assumes that the alphabet order is ACGT 
	 * @param observations
	 */
	public static void writeFasta(TreeObservations observations) {
		StringBuilder result = new StringBuilder();
		for (TreeNode node : observations.getObservedTreeNodes()) {
			result.append(">" + node + System.lineSeparator());
			double[][] s = (double[][]) observations.get(node);
			// TODO: get the alphabet order from the model
			char[] alphabet = {'A', 'C', 'G', 'T'};
			char charAtSite = 'N';
			for (int j = 0; j < s.length; j++) {
				for (int k = 0; k < s[j].length; k++) {
					if (s[j][k] != 0.0) {
						charAtSite = alphabet[k];
					}
				}
				result.append(charAtSite);
			}
			result.append(System.lineSeparator());
		}
		
	System.out.println("-----------------");
	System.out.println(result.toString());
	}
}
