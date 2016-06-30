package conifer;

import java.io.File;
import java.io.PrintWriter;

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



public class ProteinModelUsingMultipleTreesEg extends MCMCRunner
{
  File inputFile1 
    = new File("/Users/crystal/Dropbox/protein/javaProteinTry/alignment1.fasta");
  
  File inputFile2 
  = new File("/Users/crystal/Dropbox/protein/javaProteinTry/alignment2.fasta");

  @DefineFactor(onObservations = true)
  public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood1 = 
    UnrootedTreeLikelihood
    .fromFastaProteinFile(inputFile1)
    .withExpFamMixture(ExpFamMixture.proteinGTR())
    .withTree(new File("/Users/crystal/Dropbox/protein/javaProteinTry/tree1.nwk"));
  
  @DefineFactor(onObservations = true)
  public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood2 = 
    UnrootedTreeLikelihood
    .fromFastaProteinFile(inputFile2)
    .withExpFamMixture(likelihood1.evolutionaryModel.rateMatrixMixture)
    .withTree(new File("/Users/crystal/Dropbox/protein/javaProteinTry/tree2.nwk"));
  
  @DefineFactor
  NonClockTreePrior<RateParameterization> treePrior1 = 
    NonClockTreePrior
    .on(likelihood1.tree);
  
  @DefineFactor
  NonClockTreePrior<RateParameterization> treePrior2 = 
    NonClockTreePrior
    .on(likelihood2.tree);

//  @DefineFactor
//  Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = 
//    Exponential
//    .on(treePrior.branchDistributionParameters.rate)
//    .withMean(10.0);
  
  @DefineFactor
  public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior =
    IIDRealVectorGenerativeFactor
    .iidNormalOn(likelihood1.evolutionaryModel.rateMatrixMixture.parameters);
  
  private final PrintWriter treeWriter = BriefIO.output(Results.getFileInResultFolder("tree.nwk"));
  
  public static void main(String [] args)
  {
    ProteinModelUsingMultipleTreesEg runner = new ProteinModelUsingMultipleTreesEg();
    runner.factory.mcmcOptions.nMCMCSweeps = 3000;
    runner.run();
    
  }

  protected void process(ProcessorContext context)
  {
    treeWriter.println(likelihood1.tree.toNewick());
    treeWriter.flush();
  }

}
