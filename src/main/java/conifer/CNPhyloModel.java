package conifer;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.annotations.DefineFactor;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import briefj.BriefIO;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;
import conifer.ctmc.cnv.CopyNumberMixture;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.CNMultiCategorySubstitutionModel;

public class CNPhyloModel implements Runnable, Processor {
	@Option(required = true, gloss = "file containing raw reads")
	public String emissionData;

	@OptionSet(name = "factory")
	public final MCMCFactory factory = new MCMCFactory();

	@Option
	public int nMCMCSweeps = 100;

	public class Model {
		File inputFile = new File(emissionData);

		@DefineFactor(onObservations = true)
		public final UnrootedTreeLikelihood<CNMultiCategorySubstitutionModel<CopyNumberMixture>> likelihood = UnrootedTreeLikelihood
				.fromCNFile(inputFile);

		@DefineFactor
		NonClockTreePrior<RateParameterization> treePrior = NonClockTreePrior.on(likelihood.tree);

		@DefineFactor
		Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = Exponential.on(
				treePrior.branchDistributionParameters.rate).withMean(10.0);

		// TODO: maybe later put all these in a vector
		// TODO: is this a good choice of prior?
		@DefineFactor
		Exponential<Exponential.MeanParameterization> priorAlpha = Exponential.on(
				likelihood.evolutionaryModel.rateMatrixMixture.parameters.alpha).withMean(10.0);

		@DefineFactor
		Exponential<Exponential.MeanParameterization> priorBeta = Exponential.on(
				likelihood.evolutionaryModel.rateMatrixMixture.parameters.beta).withMean(10.0);

		@DefineFactor
		Exponential<Exponential.MeanParameterization> priorGamma = Exponential.on(
				likelihood.evolutionaryModel.rateMatrixMixture.parameters.gamma).withMean(10.0);
	}

	public Model model;

	private PrintWriter treeWriter;

	@Override
	public void run() {
		treeWriter = BriefIO.output(Results.getFileInResultFolder("cntrees.newick"));

		// init the model
		factory.addProcessor(this);
		model = new Model();
		MCMCAlgorithm mcmc = factory.build(model, false);
		mcmc.options.CODA = false;
		// System.out.println(mcmc.model.toString()); // TOO LONG!

		// run
		mcmc.run();

		// compute the tree
		MajorityRuleTree.buildAndWriteConsensusTree(Results.getFileInResultFolder("cntrees.newick"),
				Results.getFileInResultFolder("cnConsensusTree.Nexus"));

	}

	public static String[] getDummyArguments() {
		List<String> arguments = Arrays.asList("-emissionData",
				"src/main/resources/conifer/sampleInput/central_data_points.tsv", "-factory.mcmc.nMCMCSweeps", "100",
				"-factory.mcmc.CODA", "false", "-factory.mcmc.burnIn", "10");
		// arguments = Arrays.asList("-help");
		return arguments.toArray(new String[0]);
	}

	public static void main(String[] args) throws ClassNotFoundException, IOException {
		args = CNPhyloModel.getDummyArguments();

		for (int i = 0; i < args.length; i++)
			System.out.println(args[i]);

		Mains.instrumentedRun(args, new CNPhyloModel());
	}

	@Override
	public void process(ProcessorContext context) {
		treeWriter.println(model.likelihood.tree.toNewick());
		treeWriter.flush();
	}

}
