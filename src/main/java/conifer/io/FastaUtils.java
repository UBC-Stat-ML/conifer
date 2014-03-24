package conifer.io;

import java.io.File;
import java.util.LinkedHashMap;

import briefj.BriefIO;

import com.google.common.collect.Maps;

import conifer.TreeNode;



public class FastaUtils
{
  public static LinkedHashMap<TreeNode, CharSequence> readFasta(File f)
  {
    LinkedHashMap<TreeNode,CharSequence> map = Maps.newLinkedHashMap();
    StringBuilder current = null;
    int len = -1;
    for (String line : BriefIO.readLines(f))
      if (line.length() > 0)
      {
        if (line.charAt(0) == '>')
        {
          if (len == -1 && current != null)
            len = current.length();
          if (current != null && current.length() != len)
            throw new RuntimeException();
          
          current = new StringBuilder();
          TreeNode currentTaxon = TreeNode.withLabel(line.substring(1));
          map.put(currentTaxon, current);
        }
        else
          current.append(line);
      }
    return map;
  }
}
