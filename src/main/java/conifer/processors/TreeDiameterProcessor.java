package conifer.processors;

import java.util.List;

import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
import blang.processing.NodeProcessor;
import blang.processing.ProcessorContext;
import blang.variables.RealValued;



public class TreeDiameterProcessor implements NodeProcessor<UnrootedTree>, RealValued
{
  private UnrootedTree tree;
  private double currentValue;

  @Override
  public void process(ProcessorContext context)
  {
    this.currentValue = tree.allTotalBranchLengthDistances().max();
  }

  @Override
  public void setReference(UnrootedTree variable)
  {
    this.tree = variable;
  }

  @Override
  public double getValue()
  {
    return currentValue;
  }

}
