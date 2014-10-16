package conifer;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.io.FileUtils;

import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.MCMCRunner;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.mcmc.Move;
import blang.processing.ProcessorContext;
import briefj.BriefIO;
import briefj.opt.Option;
import briefj.run.Results;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.RateMtxNames;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;



public class SimplePhyloModel extends MCMCRunner
{
  File inputFile 
  = new File("src/main/resources/conifer/sampleInput/FES_4.fasta");

  @Option()
  public static RateMtxNames selectedRateMtx = RateMtxNames.KIMURA1980;

  @DefineFactor(onObservations = true)
  public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = 
  UnrootedTreeLikelihood
  .fromFastaFile(inputFile, selectedRateMtx)
  .withExpFamMixture(ExpFamMixture.rateMtxModel(selectedRateMtx))
  .withTree(new File("src/main/resources/conifer/sampleInput/FES.ape.4.nwk"));

  @DefineFactor
  NonClockTreePrior<RateParameterization> treePrior = 
  NonClockTreePrior
  .on(likelihood.tree);

  @DefineFactor
  Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = 
  Exponential.on(treePrior.branchDistributionParameters.rate)
  .withMean(10.0);

  @DefineFactor
  public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior =
  IIDRealVectorGenerativeFactor
  .iidNormalOn(likelihood.evolutionaryModel.rateMatrixMixture.parameters);

  private final PrintWriter treeWriter = BriefIO.output(Results.getFileInResultFolder("FES.trees.newick"));

  private final PrintWriter detailWriter = BriefIO.output(Results.getFileInResultFolder("experiment.details.txt"));
	public static void main(String [] args) throws ClassNotFoundException, IOException
	{
		SimplePhyloModel runner = new SimplePhyloModel();
		runner.factory.mcmcOptions.nMCMCSweeps = 10000;
		runner.factory.mcmcOptions.burnIn = (int) Math.round(.1 * runner.factory.mcmcOptions.nMCMCSweeps);
		
		// run
		runner.run();

  public static void main(String [] args) throws ClassNotFoundException, IOException
  {

    //		//		1.exclude: {{SingleNNI}, {SingleBranch}, {SPR}}
    //		//		
    //		//		2. exclude {{SingleNNI, SingleBranch}}
    //		//
    //		//		3. exclude {{SingleNNI, SPR}}
    //		//		
    //		//		4. exclude {{SingleBranch, SPIR}}
    //		//	
    //		//		5. exclude {{SingleNNI, SingleBranch, SPR}}
    //		//		
    //		//		6. exclude {{SPR, AllBranch}}
    @SuppressWarnings("unchecked")
    List<List<String>> experiments = Arrays.asList(Arrays.asList("SingleNNI"), 
        Arrays.asList("SingleBranchScaling"), 
        Arrays.asList("SPRMove"),
        Arrays.asList("SingleNNI", "SingleBranchScaling"),
        Arrays.asList("SingleNNI", "SPRMove"),
        Arrays.asList("SingleBranchScaling", "SPRMove"),
        Arrays.asList("SingleNNI", "SingleBranchScaling", "SPRMove"),
        Arrays.asList("SPRMove", "AllBranchesScaling")
        );
    int index = 0;
    for (List<String> experiment : experiments) {
      index++;
      long startTime = System.currentTimeMillis();

      SimplePhyloModel runner = new SimplePhyloModel();
      runner.factory.mcmcOptions.nMCMCSweeps = 1000;
      runner.factory.mcmcOptions.burnIn = (int) Math.round(.1 * runner.factory.mcmcOptions.nMCMCSweeps);
		// compute the tree
		MajorityRuleTree.buildAndWriteConsensusTree(
				Results.getFileInResultFolder("FES.trees.newick"),
				Results.getFileInResultFolder("FESConsensusTree.Nexus"));
	}

      String excluding = "Excluding " + experiment.toString();
      for (String moveClassName : experiment) {
        @SuppressWarnings("rawtypes")
        Class moveClass = Class.forName("conifer.moves." + moveClassName);
        @SuppressWarnings("unchecked")
        Class<? extends Move> castedMoveClass = (Class<? extends Move>) moveClass;
        runner.factory.excludeNodeMove(castedMoveClass);
      }
	protected void process(ProcessorContext context)
	{
		treeWriter.println(likelihood.tree.toNewick());
		treeWriter.flush();
	}

      // log experiment information
      runner.logToFile("Experiment Title:" + excluding);
      runner.logToFile("Current Date:" + RunFacility.getCurrentDateString());
      runner.logToFile("Over File:" + runner.getInputFile());
      runner.logToFile("");
      runner.logToFile("MCMC burnIn:\n" + runner.factory.mcmcOptions.burnIn);
      runner.logToFile("MCMC nMCMCSweeps:\n" + runner.factory.mcmcOptions.nMCMCSweeps);
      runner.logToFile("MCMC thinningPeriod:\n" + runner.factory.mcmcOptions.thinningPeriod);
      runner.logToFile("");
      runner.logToFile("Model:\n" + runner.buildMCMCAlgorithm().model);

      // run
      runner.run();

      //			// log finish duration
      runner.logToFile("Total time in minutes: " + ((System.currentTimeMillis() - startTime)/60000.0));
      //
      //			// compute the tree
      MajorityRuleTree.buildAndWriteConsensusTree(
          Results.getFileInResultFolder("FES.trees.newick"),
          Results.getFileInResultFolder("FESConsensusTree.Nexus"));

      // copy the results to another folder 
      File newDirectory = new File(Results.getResultFolder().getParent() + "/experiment." + Results.getResultFolder().getName() + "." + index + "." + System.currentTimeMillis());
      newDirectory.mkdir();
      FileUtils.copyDirectory(Results.getResultFolder(), newDirectory);
    }
  }

  protected void process(ProcessorContext context)
  {
    treeWriter.println(likelihood.tree.toNewick());
    treeWriter.flush();
  }

  public void logToFile(String someline) {
    this.detailWriter.println(someline);
    this.detailWriter.flush();
  }

  public String getInputFile() {
    return inputFile.getAbsolutePath();
  }
}
