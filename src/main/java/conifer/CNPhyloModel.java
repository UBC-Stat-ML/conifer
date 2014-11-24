package conifer;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

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
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;

public class CNPhyloModel implements Runnable, Processor {
	@Option(required = true, gloss = "file containing raw reads")
	public String emissionData;

	@OptionSet(name = "factory")
	public final MCMCFactory factory = new MCMCFactory();

	public class Model {
		File inputFile = new File(emissionData);

		@DefineFactor
		public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<CopyNumberMixture>> likelihood = UnrootedTreeLikelihood
				.fromCNFile(inputFile);
			
//		@DefineFactor
//		NonClockTreePrior<RateParameterization> treePrior = NonClockTreePrior.on(likelihood.tree);
//
//		@DefineFactor
//		Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = Exponential.on(
//				treePrior.branchDistributionParameters.rate).withMean(10.0);

		// emission process beta binomial precision parameters
//		@DefineFactor
//		Exponential<RateParameterization> priorGamma = Exponential.on(
//		        ((CopyNumberTreeObservation) likelihood.observations).betaBinomialprecision);
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
		File graph = Results.getFileInResultFolder("probability-graph.dot");	
		mcmc.model.printGraph(graph);
		mcmc.options.CODA = true;
		System.out.println(mcmc.model.toString());

		// run
		mcmc.run();

		// compute the tree
		MajorityRuleTree.buildAndWriteConsensusTree(Results.getFileInResultFolder("cntrees.newick"),
				Results.getFileInResultFolder("cnConsensusTree.Nexus"));

	}

	public static void main(String[] args) throws ClassNotFoundException, IOException {
		Mains.instrumentedRun(args, new CNPhyloModel());
	}

	@Override
	public void process(ProcessorContext context) {
		treeWriter.println(model.likelihood.tree.toNewick());
		treeWriter.flush();
	}

}
