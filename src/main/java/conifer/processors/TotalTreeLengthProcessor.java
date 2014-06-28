package conifer.processors;

import conifer.UnrootedTree;
import conifer.UnrootedTreeUtils;
import blang.processing.NodeProcessor;
import blang.processing.ProcessorContext;
import blang.variables.RealValued;


/**
 * Computes the sum of the branch lengths of a tree.
 * Used for testing.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class TotalTreeLengthProcessor implements NodeProcessor<UnrootedTree>, RealValued
{
  private UnrootedTree tree;
  private double currentValue;

  @Override
  public void process(ProcessorContext context)
  {
    currentValue = UnrootedTreeUtils.totalTreeLength(tree);
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
