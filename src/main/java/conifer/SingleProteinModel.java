package conifer;

import java.io.File;
import java.io.PrintWriter;
import java.util.Random;

//import bayonet.distributions.Exponential;
//import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.ForwardSampler;
import blang.MCMCRunner;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.processing.ProcessorContext;
import briefj.BriefIO;
import briefj.opt.Option;
import briefj.opt.OptionsParser;
import briefj.tomove.Results;
import conifer.ctmc.expfam.ExpFamMixture;
//import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;



public class SingleProteinModel extends MCMCRunner
{
  
 // private String myPath;
  public static class ProteinOptions
  {
   @Option(gloss="File of provided alignment")
   public static File inputFile;
   @Option(gloss="File of the tree topology")
   public static File treeFile;
  //File inputFile = new File("/Users/crystal/Dropbox/protein/data/titin/alignmentContact.txt");
  }
  @DefineFactor(onObservations = true)
  public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood1 = 
    UnrootedTreeLikelihood
    .fromFastaProteinFile(ProteinOptions.inputFile)
    .withExpFamMixture(ExpFamMixture.proteinGTR())
    .withTree(ProteinOptions.treeFile);
   // .withTree(new File("/Users/crystal/Dropbox/protein/data/titin/tree.nwk.txt"));
  
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
  
  @DefineFactor
  public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior =
    IIDRealVectorGenerativeFactor
    .iidNormalOn(likelihood1.evolutionaryModel.rateMatrixMixture.parameters);
  
 // private final PrintWriter treeWriter = BriefIO.output(Results.getFileInResultFolder("tree.nwk"));
  
  public static void main(String [] args)
  {
    
    OptionsParser optionParser = new OptionsParser();
    optionParser.register(ProteinOptions.class);
    optionParser.parse(args);
    
    
    
    SingleProteinModel runner = new SingleProteinModel();
    
    ProbabilityModel model = ProbabilityModel.parse(runner);
    ForwardSampler sampler = new ForwardSampler(model);
    Random rand = new Random(1);
    sampler.simulate(rand);
    
    
    //runner.setMyPath("");
    runner.factory.mcmcOptions.nMCMCSweeps = 100;
    runner.run();
  }

 /* public String getMyPath()
  {
    return myPath;
  }

  public void setMyPath(String myPath)
  {
    this.myPath = myPath;
  }
*/
//  protected void process(ProcessorContext context)
  //{
  //  treeWriter.println(likelihood1.tree.toNewick());
   // treeWriter.flush();
  //}

}
