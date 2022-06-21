package conifer;

import org.junit.Assert;
import org.junit.Test;

import briefj.BriefIO;



public class CladeTest
{
  @Test
  public void test()
  {
//    if (true)
//      throw new RuntimeException();
    UnrootedTree tree = UnrootedTreeUtils.fromNewickString(BriefIO.resourceToString("/conifer/smallTree.newick", CladeTest.class));
    Assert.assertEquals(UnrootedTreeUtils.bipartitions(tree).size(), 5);
  }
}
