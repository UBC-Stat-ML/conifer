package conifer.io;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import org.apache.commons.io.FileUtils;

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
	 * @throws IOException 
	 */
	public static String writeFasta(TreeObservations observations, File outFile) throws IOException 
	{
		// get the alphabet map
		Map<String,String> a2s = PhylogeneticObservationFactory.nucleotidesFactory().getIndicator2ChunkMap();

		StringBuilder result = new StringBuilder();
		for (TreeNode node : observations.getObservedTreeNodes()) {
			result.append(">" + node + "\n");
			double[][] s = (double[][]) observations.get(node);
			char charAtSite = 'N';
			for (int j = 0; j < s.length; j++) {
				charAtSite = a2s.get(Arrays.toString(s[j])).charAt(0);
				result.append(charAtSite);
			}
			result.append("\n");
		}
		
		FileUtils.writeStringToFile(outFile, result.toString());

		return result.toString();
	}


}
