package conifer.processors;

import blang.processing.NodeProcessor;
import blang.processing.ProcessorContext;
import blang.variables.RealValued;
import conifer.UnrootedTree;


/**
 * The tree diameter is the max distance between two leaves of a tree.
 * Used for testing.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
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
