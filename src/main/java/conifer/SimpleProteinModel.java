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
import blang.processing.Processor;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Results;
import briefj.run.Mains;
import briefj.BriefIO;
import briefj.OutputManager;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;

public class SimpleProteinModel implements Runnable, Processor
{

  @Option(gloss="File of provided alignment")
  public static  File inputFile1 
  = new File("/Users/crystal/Dropbox/protein/javaProteinTry/alignment1.fasta");

  @Option(gloss="File of provided alignment")
  public static File inputFile2 
  = new File("/Users/crystal/Dropbox/protein/javaProteinTry/alignment2.fasta");

  @Option(gloss="File of the tree topology")
  public static File treeFile1 = new File("/Users/crystal/Dropbox/protein/javaProteinTry/tree1.nwk");

  @Option(gloss="File of the tree topology")
  public static File treeFile2 = new File("/Users/crystal/Dropbox/protein/javaProteinTry/tree2.nwk");



  @Option()
  public static String selectedRateMtx="polaritySize()";

  @OptionSet(name = "factory")
  public final MCMCFactory factory = new MCMCFactory();


  public class Model
  {
    @DefineFactor(onObservations = true)
    public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood1 = 
    UnrootedTreeLikelihood
    .fromFastaFile(inputFile1, selectedRateMtx)
    .withExpFamMixture(ExpFamMixture.rateMtxModel(selectedRateMtx))
    .withTree(treeFile1);

    @DefineFactor(onObservations = true)
    public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood2 = 
    UnrootedTreeLikelihood
    .fromFastaFile(inputFile2, selectedRateMtx)
    .withExpFamMixture(likelihood1.evolutionaryModel.rateMatrixMixture)
    .withTree(treeFile2);

    @DefineFactor
    NonClockTreePrior<RateParameterization> treePrior1 = 
    NonClockTreePrior
    .on(likelihood1.tree);

    @DefineFactor
    NonClockTreePrior<RateParameterization> treePrior2 = 
    NonClockTreePrior
    .on(likelihood2.tree);

    // @DefineFactor
    // Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = 
    //   Exponential
    //   .on(treePrior.branchDistributionParameters.rate)
    //   .withMean(10.0);


    @DefineFactor
    public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior =
    IIDRealVectorGenerativeFactor
    .iidNormalOn(likelihood1.evolutionaryModel.rateMatrixMixture.parameters);
  }

  public Model model;
  private final PrintWriter treeWriter = BriefIO.output(Results.getFileInResultFolder("tree.nwk"));

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
    Mains.instrumentedRun(args, new SimpleProteinModel());
  }



  public void process(ProcessorContext context)
  {
    treeWriter.println(model.likelihood1.tree.toNewick());
    treeWriter.flush();
  }


}













  
  
  
  


  

  
 


