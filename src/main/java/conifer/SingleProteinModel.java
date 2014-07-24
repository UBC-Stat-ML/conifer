package conifer;

import java.io.File;
import java.util.Random;

import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.ForwardSampler;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.MCMCRunner;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.processing.Processor;
import briefj.opt.Option;
import briefj.opt.OptionsParser;
import briefj.opt.OptionSet;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.RateMtxNames;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import briefj.run.Mains;
import blang.processing.ProcessorContext;
import briefj.OutputManager;



public class SingleProteinModel implements Runnable, Processor
{
   @Option(gloss="File of provided alignment")
   public File inputFile;

   @Option(gloss="File of the tree topology")
   public File treeFile;
   
   @Option()
   public static RateMtxNames selectedRateMtx = RateMtxNames.POLARITYSIZE;

   @OptionSet(name = "factory")
   public final MCMCFactory factory = new MCMCFactory();
  
  
  public class Model
  {
    @DefineFactor(onObservations = true)
    public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood1 = 
    UnrootedTreeLikelihood
    .fromFastaFile(inputFile, selectedRateMtx)
    .withExpFamMixture(ExpFamMixture.rateMtxModel(selectedRateMtx))
    .withTree(treeFile);

    @DefineFactor
    public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior =
    IIDRealVectorGenerativeFactor
    .iidNormalOn(likelihood1.evolutionaryModel.rateMatrixMixture.parameters);
  }
  
  
// The topology of the tree is fixed so that I don't put a prior on the tree topology
// @DefineFactor
// NonClockTreePrior<RateParameterization> treePrior1 = 
//  NonClockTreePrior
// .on(likelihood1.tree);
  
  
//  @DefineFactor
//  Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = 
//    Exponential
//    .on(treePrior.branchDistributionParameters.rate)
//    .withMean(10.0);
  
// private final PrintWriter treeWriter = BriefIO.output(Results.getFileInResultFolder("tree.nwk"));
 
// Note: only instantiate this in run to avoid problems with command line argument parsing
  public Model model;

  public void run()
  {
    factory.addProcessor(this);
    model = new Model();
    MCMCAlgorithm mcmc = factory.build(model, false);
    mcmc.options.nMCMCSweeps=10000; 
    mcmc.run();
  }

 public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new SingleProteinModel());
  }

  @Override
  public void process(ProcessorContext context)
  {
    
  }

}
