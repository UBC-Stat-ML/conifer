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
import briefj.opt.Option;
import briefj.run.Results;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.RateMtxNames;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.moves.AllBranchesScaling;

public class HierarchicalPhyloModel extends MCMCRunner {
	private boolean useOldMove = false;

	File inputFileUTY = new File(
			"src/main/resources/conifer/sampleInput/UTY_4.fasta");

	File inputFileFES
	= new File(
			"src/main/resources/conifer/sampleInput/FES_4.fasta");

	@Option(gloss="Indicator of we normalize the rate matrix if it is set to true")
	public boolean isNormalized = true;
	
	@Option()
  public static RateMtxNames selectedRateMtx=RateMtxNames.KIMURA1980;

	@DefineFactor(onObservations = true)
	public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihoodUTY = UnrootedTreeLikelihood
	.fromFastaFile(inputFileUTY, selectedRateMtx)
	.withExpFamMixture(ExpFamMixture.rateMtxModel(selectedRateMtx, isNormalized))
	.withTree(new File("src/main/resources/conifer/sampleInput/UTY.ape.4.nwk"));

	@DefineFactor(onObservations = true)
	public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihoodFES = UnrootedTreeLikelihood
	.fromFastaFile(inputFileFES, selectedRateMtx)
	.withExpFamMixture(likelihoodUTY.evolutionaryModel.rateMatrixMixture)
	.withTree(new File("src/main/resources/conifer/sampleInput/FES.ape.4.nwk"));

	@DefineFactor
	NonClockTreePrior<RateParameterization> treePriorUTY = NonClockTreePrior
	.on(likelihoodUTY.tree);

	@DefineFactor
	NonClockTreePrior<RateParameterization> treePriorFES = NonClockTreePrior
	.on(likelihoodFES.tree);

	// define hyper priors
	@DefineFactor
	Exponential<Exponential.RateParameterization> branchLengthHyperPriorUTY = Exponential
	.on(treePriorUTY.branchDistributionParameters.rate).withRate(0.1);

	// define hyper-hyper prior
	@DefineFactor
	Exponential<Exponential.RateParameterization> branchLengthHyperHyperPriorUTY = Exponential
	.on(branchLengthHyperPriorUTY.parameters.rate)
	.withRate(0.1);

	// share the prior between trees
	@DefineFactor
	Exponential<Exponential.RateParameterization> branchLengthHyperPriorFES = Exponential
	.on(treePriorFES.branchDistributionParameters.rate)
	.withRate(branchLengthHyperPriorUTY.parameters.rate);

	@DefineFactor
	public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior = IIDRealVectorGenerativeFactor
	.iidNormalOn(likelihoodUTY.evolutionaryModel.rateMatrixMixture.parameters);

	private final PrintWriter treeWriterFES = BriefIO.output(Results
			.getFileInResultFolder("FES.trees.newick"));
	private final PrintWriter treeWriterUTY = BriefIO.output(Results
			.getFileInResultFolder("UTY.trees.newick"));

	private final PrintWriter detailWriter = BriefIO.output(Results.getFileInResultFolder("experiment.details.txt"));


	public boolean isUseOldMove() {
		return useOldMove;
	}

	public void setUseOldMove(boolean useOldMove) {
		this.useOldMove = useOldMove;
	}

	public static void main(String[] args) {

		HierarchicalPhyloModel runner = new HierarchicalPhyloModel();
		runner.factory.mcmcOptions.nMCMCSweeps = 10000;
		runner.factory.mcmcOptions.burnIn = (int) Math.round(.1 * runner.factory.mcmcOptions.nMCMCSweeps);

		String excluding = "All moves, hierarchical, sharing rate matrix and hyperhyper prior. 8 Taxa";

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

		
		runner.factory.excludeNodeMove(AllBranchesScaling.class);
		// run
		runner.run();

		MajorityRuleTree.buildAndWriteConsensusTree(
				Results.getFileInResultFolder("UTY.trees.newick"),
				Results.getFileInResultFolder("UTYConsensusTree.Nexus"));

		MajorityRuleTree.buildAndWriteConsensusTree(
				Results.getFileInResultFolder("FES.trees.newick"),
				Results.getFileInResultFolder("FESConsensusTree.Nexus"));

	}

	protected void process(ProcessorContext context) {
		treeWriterFES.println(likelihoodFES.tree.toNewick());
		treeWriterFES.flush();

		treeWriterUTY.println(likelihoodUTY.tree.toNewick());
		treeWriterUTY.flush();

		// hyperPriorWriterFES.println(branchLengthHyperHyperPriorUTY.parameters.rate);
		// hyperPriorWriterFES.println("lambda_0="+branchLengthHyperPriorUTY.parameters.rate);
		// hyperPriorWriterFES.println("lambda_0="+branchLengthHyperPriorFES.parameters.rate);
	}

	public void logToFile(String someline) {
		this.detailWriter.println(someline);
		this.detailWriter.flush();
	}

	public String getInputFile() {
		return inputFileFES.getAbsolutePath() + "\n" + inputFileUTY.getAbsolutePath();
	}

}
