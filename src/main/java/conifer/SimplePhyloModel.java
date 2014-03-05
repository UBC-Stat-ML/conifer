package conifer;

import java.io.File;
import java.util.List;

import com.google.common.collect.Lists;

import conifer.ctmc.JukeCantorRateMatrix;
import conifer.factors.NonClockTreePrior;
import conifer.factors.TreeLikelihood;
import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import bayonet.rplot.PlotHistogram;
import blang.MCMCRunner;
import blang.annotations.DefineFactor;
import blang.processing.ProcessorContext;
import briefj.tomove.Results;



public class SimplePhyloModel extends MCMCRunner
{
  File inputFile = new File("primates.data");
  
  @DefineFactor 
  TreeLikelihood<JukeCantorRateMatrix> treeLikelihood = TreeLikelihood.fromObservations(inputFile);
  
  @DefineFactor
  NonClockTreePrior<Exponential<RateParameterization>> treePrior = NonClockTreePrior.on(treeLikelihood.tree);

  @DefineFactor
  Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = Exponential.on(treePrior.branchDistribution.parameters.rate).withMean(10.0);
  
  
  
  
  public static void main(String [] args)
  {
    SimplePhyloModel runner = new SimplePhyloModel();
    runner.factory.mcmcOptions.nMCMCSweeps = 10000;
    runner.run();
    File plotFile = Results.getFileInResultFolder("density.pdf");
    PlotHistogram.from(runner.samples).toPDF(plotFile);
  }


  protected void process(ProcessorContext context)
  {
    samples.add(treePrior.branchDistribution.parameters.rate.getValue());
  }
  
  List<Double> samples = Lists.newArrayList();
}
