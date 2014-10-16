package conifer.processors;

import java.util.List;

import blang.MCMCAlgorithm;
import blang.processing.NodeProcessor;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import blang.variables.RealValued;
import blang.variables.RealVariableProcessor;

import com.google.common.collect.Lists;

import conifer.UnrootedTree;
import conifer.UnrootedTreeUtils;
import conifer.UnrootedTreeUtils.TreeMetric;


public class TreeDistanceProcessor implements NodeProcessor<UnrootedTree>
{
  public static void setReferenceTree(MCMCAlgorithm mcmc, UnrootedTree sampledTree, UnrootedTree referenceTree)
  {
    for (Processor p : mcmc.processors)
      if (p instanceof TreeDistanceProcessor)
      {
        TreeDistanceProcessor current = (TreeDistanceProcessor) p;
        if (current.sampledTree == sampledTree)
          current.setReferenceTree(referenceTree);
      }
  }
  
  private UnrootedTree referenceTree, sampledTree;
  private List<RealVariableProcessor> subProcessors = null;

  private void setReferenceTree(UnrootedTree referenceTree)
  {
    this.referenceTree = new UnrootedTree(referenceTree);
  }
  
  @Override
  public void process(ProcessorContext context)
  {
    if (!ensureInitialized(context))
      return; // if reference not set, do nothing
    
    for (RealVariableProcessor subprocessor : subProcessors)
      subprocessor.process(context);
  }


  private boolean ensureInitialized(ProcessorContext context)
  {
    if (referenceTree == null)
      return false;
    if (subProcessors == null)
    {
      subProcessors = Lists.newArrayList();
      String treeVariableName = context.getModel().getName(sampledTree);
      for (TreeMetric metric : TreeMetric.values())
        subProcessors.add(new RealVariableProcessor(treeVariableName + "-" + metric, calculator(metric)));
    }
    return true;
  }

  private RealValued calculator(final TreeMetric metric)
  {
    return new RealValued() {
      
      @Override
      public double getValue()
      {
        return UnrootedTreeUtils.computeTreeMetrics(referenceTree, sampledTree).get(metric);
      }
    };
  }

  @Override
  public void setReference(UnrootedTree variable)
  {
    this.sampledTree = variable; 
  }
}
