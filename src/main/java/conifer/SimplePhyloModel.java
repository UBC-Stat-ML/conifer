package conifer;

import java.io.File;
import java.io.PrintWriter;

import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.MCMCRunner;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.processing.ProcessorContext;
import briefj.BriefIO;
import briefj.run.Results;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;



public class SimplePhyloModel extends MCMCRunner
{
  File inputFile 
//    = new File("primates.fasta");
    = new File("/Users/bouchard/Documents/data/utcs/23S.E/R0/cleaned.alignment.fasta");
  
  @DefineFactor(onObservations = true)
  public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = 
    UnrootedTreeLikelihood
    .fromFastaFile(inputFile)
    .withExpFamMixture(ExpFamMixture.kimura1980())
    .withTree(new File("/Users/bouchard/Documents/data/utcs/23S.E.raxml.nwk"));
  
  @DefineFactor
  NonClockTreePrior<RateParameterization> treePrior = 
    NonClockTreePrior
    .on(likelihood.tree);

  @DefineFactor
  Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = 
    Exponential
    .on(treePrior.branchDistributionParameters.rate)
    .withMean(10.0);
  
  @DefineFactor
  public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior =
    IIDRealVectorGenerativeFactor
    .iidNormalOn(likelihood.evolutionaryModel.rateMatrixMixture.parameters);
  
  private final PrintWriter treeWriter = BriefIO.output(Results.getFileInResultFolder("trees.newick"));
  
  public static void main(String [] args)
  {
    SimplePhyloModel runner = new SimplePhyloModel();
    runner.factory.mcmcOptions.nMCMCSweeps = 10000000;
    runner.run();
    
  }

  protected void process(ProcessorContext context)
  {
    treeWriter.println(likelihood.tree.toNewick());
    treeWriter.flush();
  }

}
