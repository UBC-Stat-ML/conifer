package conifer.io.newick;


import org.junit.Assert;
import org.junit.Test;

import conifer.UnrootedTree;
import conifer.UnrootedTreeUtils;



public class TestParser
{
  @Test
  public void test()
  {
    String treeStr1 = "((A:7.0,B:6.0):5.0,(C:4.0,D:3.0)I:2.0,E:1.0);";
    UnrootedTree tree1 = UnrootedTree.fromNewickString(treeStr1);
    Assert.assertEquals(treeStr1, tree1.toString());
    
    String treeStr2 = "((C:4.0,D:3.0)I:2.0,(A:7.0,B:6.0):5.0,E:1.0);";
    UnrootedTree tree2 = UnrootedTree.fromNewickString(treeStr2);
    Assert.assertEquals(UnrootedTreeUtils.computeTreeMetrics(tree1, tree2).get(UnrootedTreeUtils.TreeMetric.l1), 0.0, 1e-10);
    
    Assert.assertEquals(28.0, UnrootedTreeUtils.totalTreeLength(tree2), 1e-10);
  }
}
