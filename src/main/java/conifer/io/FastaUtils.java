package conifer.io;

import java.io.File;
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
}
