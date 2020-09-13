package conifer;

import java.io.File;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.rits.cloning.Immutable;

import blang.inits.ConstructorArg;
import blang.inits.DesignatedConstructor;
import blang.inits.GlobalArg;
import blang.runtime.Observations;
import blang.runtime.internals.objectgraph.SkipDependency;
import briefj.BriefCollections;
import conifer.io.FastaUtils;
import conifer.io.PhylogeneticObservationFactory;
import conifer.io.TreeObservations;

/**
 * By Fixed, we mean that the alignment is not altered. 
 */
public class SequenceAlignment implements TreeObservations  
{
  @SkipDependency(isMutable = true)
  private final LinkedHashMap<TreeNode, double[][]> data = Maps.newLinkedHashMap();
  
  private final int nSites;
  
  @SkipDependency(isMutable = false)
  public final PhylogeneticObservationFactory factory;
  
  public SequenceAlignment(PhylogeneticObservationFactory factory, int nSites) 
  {
    this.factory = factory;
    this.nSites = nSites;
  }
  
  @DesignatedConstructor
  public static SequenceAlignment loadObservedData( 
      @ConstructorArg(value = "file",     description = "Path to FASTA alignment")             File                           file,
      @ConstructorArg(value = "encoding", description = "Encoding used to read the alignment") PhylogeneticObservationFactory factory,
      @GlobalArg Observations observations)
  {
    LinkedHashMap<TreeNode, String> readFasta = FastaUtils.readFasta(file); 
    SequenceAlignment result = new SequenceAlignment(factory, BriefCollections.pick(readFasta.values()).length());
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
      throw new RuntimeException("Make sure the data is aligned. Number of sites seems to differ for different leaves: " + cast.length + " vs " + nSites);
    
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
