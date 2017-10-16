package conifer.io;

import java.io.File;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import blang.inits.ConstructorArg;
import blang.inits.DesignatedConstructor;
import blang.inits.GlobalArg;
import blang.runtime.Observations;
import briefj.BriefCollections;
import conifer.TreeNode;

/**
 * By Fixed, we mean that the alignment is not altered. 
 */
public class FixedTreeObservations implements TreeObservations  
{
  private final LinkedHashMap<TreeNode, double[][]> data = Maps.newLinkedHashMap();
  private final int nSites;
  
  public final PhylogeneticObservationFactory factory;
  
  public FixedTreeObservations(PhylogeneticObservationFactory factory, int nSites) 
  {
    this.factory = factory;
    this.nSites = nSites;
  }
  
  @DesignatedConstructor
  public static FixedTreeObservations loadObservedData( 
      @ConstructorArg(value = "file",     description = "Path to FASTA alignment")             File                           file,
      @ConstructorArg(value = "encoding", description = "Encoding used to read the alignment") PhylogeneticObservationFactory factory,
      @GlobalArg Observations observations)
  {
    LinkedHashMap<TreeNode, String> readFasta = FastaUtils.readFasta(file); 
    FixedTreeObservations result = new FixedTreeObservations(factory, BriefCollections.pick(readFasta.values()).length());
    for (TreeNode treeNode : readFasta.keySet())
    {
        String rawString = readFasta.get(treeNode); 
        double [][] indicators = factory.site2CharacterIndicators(rawString);
        result.set(treeNode, indicators);
    }
    // Note: even sites with uncertainty can be marked as observed since their uncertainty will be marginalized analytically.
    observations.markAsObserved(result);
    return result;
  }

  @Override
  public List<TreeNode> getObservedTreeNodes()
  {
    List<TreeNode> result = Lists.newArrayList();
    for (TreeNode key : data.keySet())
      result.add(key);
    return result;
  }

  @Override
  public double[][] get(TreeNode leaf)
  {
    return data.get(leaf);
  }

  @Override
  public void set(TreeNode leaf, Object point)
  {
    double[][] cast = (double[][])point;
    if (cast.length != nSites)
      throw new RuntimeException();
    data.put(leaf, cast);
  }

  @Override
  public void clear()
  {
    data.clear();
  }

  @Override
  public String toString()
  {
    StringBuilder result = new StringBuilder();
    for (TreeNode node : data.keySet())
      result.append(node.toString() + " : " + Arrays.deepToString(data.get(node)) + " ");
    return result.toString();
  }

  @Override
  public int nSites()
  {
    return nSites;
  }
}
