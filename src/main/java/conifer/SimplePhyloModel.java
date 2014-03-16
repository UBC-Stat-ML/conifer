package conifer;

import java.io.File;
import java.util.List;

import tutorialj.Tutorial;
import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import bayonet.rplot.PlotHistogram;
import blang.MCMCRunner;
import blang.annotations.DefineFactor;
import blang.processing.ProcessorContext;
import briefj.tomove.Results;

import com.google.common.collect.Lists;

import conifer.ctmc.JukeCantorRateMatrix;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;



public class SimplePhyloModel extends MCMCRunner
{
  File inputFile = new File("primates.data");
  
  @DefineFactor(onObservations = true)
  UnrootedTreeLikelihood<JukeCantorRateMatrix> treeLikelihood = UnrootedTreeLikelihood.fromObservations(inputFile);
  
  @DefineFactor
  NonClockTreePrior<RateParameterization> treePrior = NonClockTreePrior.on(treeLikelihood.tree);

  @DefineFactor
  Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = Exponential.on(treePrior.branchDistributionParameters.rate).withMean(10.0);
  
  
  
  
  public static void main(String [] args)
  {
    SimplePhyloModel runner = new SimplePhyloModel();
    runner.factory.mcmcOptions.nMCMCSweeps = 10000;
    runner.run();
    File plotFile = Results.getFileInResultFolder("density.pdf");
    PlotHistogram.from(runner.samples).toPDF(plotFile);
  }

  /**
   * Modify the ``process()`` function below to save the output trees to a file.
   * Use ``UnrootedTreeUtils.toNewick()`` to save the files in the 
   * [newick format](http://en.wikipedia.org/wiki/Newick_format), a standard 
   * way of saving tree.
   */
  @Tutorial(showSource = true, showLink = true, showSignature = true)
  protected void process(ProcessorContext context)
  {
    samples.add(treePrior.branchDistributionParameters.rate.getValue());
  }
  
  List<Double> samples = Lists.newArrayList();
  

}
