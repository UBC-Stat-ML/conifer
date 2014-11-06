package conifer;

import java.io.File;
import java.io.PrintWriter;

import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import briefj.BriefIO;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;


/**
 * 
 * Build a basic phylo model
 * test the resampling of the rate matrix parameterizations
 * 
 * @author Sean Jewell (jewellsean@gmail.com)
 * CTMC
 */
public class BasicPhyloModel implements Runnable, Processor 
{
  @Option
  File inputFile 
  = new File("src/main/resources/conifer/sampleInput/FES_4.fasta");
  
  @OptionSet(name = "factory")
  public final MCMCFactory factory = new MCMCFactory();
  
  /**
   * 
   * Simple phylo model
   * 
   */
 
  public class Model
  {
    @DefineFactor(onObservations = true)
    public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = 
    UnrootedTreeLikelihood
    .fromFastaFile(inputFile)
    .withExpFamMixture(ExpFamMixture.kimura1980());

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
    public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior = IIDRealVectorGenerativeFactor
    .iidNormalOn(likelihood.evolutionaryModel.rateMatrixMixture.parameters);
  }
  
  public Model model;
  
  private final PrintWriter treeWriter =
      BriefIO.output(Results.getFileInResultFolder("FES.trees.newick"));
  
  @Override
  public void run()
  {
   model = new Model(); 
   MCMCAlgorithm mcmc = factory.build(model, false); // Q: why do we need to clone this? 
   System.out.println(mcmc.model);
   mcmc.run(); 
  }

  public static void main(String args[])
  {
     Mains.instrumentedRun(args, new BasicPhyloModel());
  }
  
  
  @Override
  public void process(ProcessorContext context){
    System.out.println(model.likelihood.tree.toNewick());
    treeWriter.flush();
  }
  
}
