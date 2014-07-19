package conifer;

import org.junit.Assert;
import org.junit.Test;

import briefj.BriefIO;



public class CladeTest
{
  @Test
  public void test()
  {
    UnrootedTree tree = UnrootedTreeUtils.fromNewickString(BriefIO.resourceToString("/conifer/smallTree.newick"));
    Assert.assertEquals(UnrootedTreeUtils.bipartitions(tree).size(), 5);
  }
}
