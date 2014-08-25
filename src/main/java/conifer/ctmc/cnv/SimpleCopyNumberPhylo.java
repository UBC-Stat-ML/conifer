//package conifer.ctmc.cnv;
//
//import java.io.File;
//
//import conifer.ctmc.expfam.ExpFamMixture;
//import conifer.factors.NonClockTreePrior;
//import conifer.factors.UnrootedTreeLikelihood;
//import conifer.models.DiscreteGammaMixture;
//import conifer.models.MultiCategorySubstitutionModel;
//import bayonet.distributions.Exponential;
//import bayonet.distributions.Exponential.RateParameterization;
//import blang.MCMCAlgorithm;
//import blang.MCMCFactory;
//import blang.annotations.DefineFactor;
//import blang.processing.Processor;
//import blang.processing.ProcessorContext;
//import briefj.opt.Option;
//import briefj.opt.OptionSet;
//import briefj.run.Mains;
//
///**
// * Implement a simple phylo model with rate matrix given by
// * sparse representation 
// * 
// * Q = mutationEpsilon * mutationQ + 
// *     increaseCNV * increaseQ + 
// *     decreaseCNV * decreaseQ 
// * 
// * @author Sean Jewell (jewellsean@gmail.com)
// *
// */
//
//public class SimpleCopyNumberPhylo implements Runnable, Processor
//{
//  @Option
//  File inputSiteData
//   = new File("...fasta"); // need some data or some way to generate it
//  
//  @OptionSet(name = "factory")
//  public final MCMCFactory factory = new MCMCFactory(); 
//  
//  
//  /**
//   * 
//   * Simple implementation 
//   * 
//   * Potential problems: 
//   *        1) Rate matrix normalizaiton working correctly 
//   *        2) Indentification problmes
//   *        3) values of increaseCNV and decreaseCNV are stored pre norm
//   *        
//   *        this seems like a problem
//   *        
//   */
//
//  public class Model
//  {
//    // add support in this likelihood
//    @DefineFactor(onObservations = true)
//    public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<CopyNumberMixture>> likelihood = 
//    UnrootedTreeLikelihood
//    .fromFastaFile(inputSiteData)
//
//    @DefineFactor
//    NonClockTreePrior<RateParameterization> treePrior = 
//    NonClockTreePrior
//    .on(likelihood.tree);
//
//    @DefineFactor
//    Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = 
//    Exponential
//    .on(treePrior.branchDistributionParameters.rate)
//    .withMean(10.0);
//    
//    @DefineFactor
//    Exponential<Exponential.MeanParameterization> increaseCNV = 
//    Exponential
//    .on(likelihood.evolutionaryModel.rateMatrixMixture.parameters.increaseCNV)
//    .withMean(10.0);
//    
//    @DefineFactor
//    Exponential<Exponential.MeanParameterization> decreaseCNV = 
//    Exponential
//    .on(likelihood.evolutionaryModel.rateMatrixMixture.parameters.decreaseCNV)
//    .withMean(10.0);
//  }
//  
//  public Model model;
//  
//  @Override
//  public void run()
//  {
//    model = new Model(); 
//    MCMCAlgorithm mcmc = factory.build(model, false);
//    System.out.println(mcmc.model);
//    mcmc.run();
//  }
//
//  public static void main(String args [])
//  {
//    Mains.instrumentedRun(args, new SimpleCopyNumberPhylo());
//  }
//  
//  @Override
//  public void process(ProcessorContext context) {}
//  
//}
