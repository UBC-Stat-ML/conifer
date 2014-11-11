package conifer;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.ForwardSampler;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import blang.variables.RealValued;
import briefj.BriefIO;
import briefj.BriefMaps;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.io.FastaUtils;
import conifer.io.TreeObservations;
import conifer.models.MultiCategorySubstitutionModel;




public class TestPhyloModel implements Runnable, Processor
{
	@Option
	public String initialTreeFilePath;
	//= new File("/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES_4.fasta");

	@Option
	public String alignmentFilePath = "";

	@Option
	public int nMCMCSweeps = 100;

	@Option
	public int burnIn = (int) Math.round(.1 * nMCMCSweeps);

	@Option
	public int thinningPeriod = 10;

	@Option
	public String phyloModel = "GTR";

	@OptionSet(name = "factory")
	public final MCMCFactory factory = new MCMCFactory();


	@Option
	public int nSites = 500;

	public class Model
	{
		List<TreeNode> leaves = Arrays.asList(
				TreeNode.withLabel("Human"), 
				TreeNode.withLabel("Chimp"),
				TreeNode.withLabel("Gorilla"), 
				TreeNode.withLabel("Orangutan"));

		@DefineFactor(onObservations = true)
		public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = 
		UnrootedTreeLikelihood//.createEmpty(nSites, leaves)
		.fromFastaFile(new File(alignmentFilePath))
		.withExpFamMixture(ExpFamMixture.kimura1980())
		.withTree(new File(initialTreeFilePath));

		@DefineFactor
		NonClockTreePrior<RateParameterization> treePrior = 
		NonClockTreePrior
		.on(likelihood.tree);

		@DefineFactor
		Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = 
		Exponential
		.on(treePrior.branchDistributionParameters.rate)
		.withMean(0.1);

		@DefineFactor
		public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior =
		IIDRealVectorGenerativeFactor
		.iidNormalOn(likelihood.evolutionaryModel.rateMatrixMixture.parameters);		
	}

	// we can build different models here and instantiate them according to the input parameter
	// in the run block.

	// Note: only instantiate this in run to avoid problems with command line argument parsing
	public Model model;

	private PrintWriter treeWriter;
	private PrintWriter detailWriter;

	@SuppressWarnings("unchecked")
	@Override
	public void run()
	{
		System.out.println("Running the newest version...");
		treeWriter = BriefIO.output(Results.getFileInResultFolder("FES.trees.newick"));
		detailWriter = BriefIO.output(Results.getFileInResultFolder("experiment.details.txt"));

		
		logToFile("Over AlignmentFile:" + getAlignmentFile());
		
		factory.addProcessor(this);
		factory.mcmcOptions.nMCMCSweeps = nMCMCSweeps;
		factory.mcmcOptions.burnIn = burnIn;
		factory.mcmcOptions.CODA = true;
		factory.mcmcOptions.thinningPeriod = thinningPeriod;

		// reference bug http://stackoverflow.com/questions/4829576/javac-error-inconvertible-types-with-generics
		
		// Fixed topology, , latent branch length
//		factory.excludeNodeMove((Class<? extends Move>) (Object)SingleNNI.class);
//		factory.excludeNodeMove(SPRMove.class);
		
		// Fixed topology, fixed branch length
		//factory.excludeNodeMove((Class<? extends Move>) (Object)SingleBranchScaling.class);
		//factory.excludeNodeMove(AllBranchesScaling.class);
		
		//factory.excludeNodeMove((Class<? extends Move>) (Object)SingleBranchScaling.class);
		//factory.excludeNodeMove(AllBranchesScaling.class);
		//factory.excludeNodeMove(PhyloHMCMove.class);
		//factory.excludeNodeMove((Class<? extends Move>) (Object)SingleNNI.class);
		//factory.excludeNodeMove((Class<? extends Move>) (Object)RealVectorMHProposal.class);
		
		model = new Model();
		MCMCAlgorithm mcmc = factory.build(model, false);
		System.out.println(mcmc.model);

		long startTime = System.currentTimeMillis();
		String excluding = "GTR-Excluding all but SPR move";

		// log experiment information
		logToFile("Experiment Title:" + excluding);
		logToFile("Current Date:" + RunFacility.getCurrentDateString());
		logToFile("Over AlignmentFile:" + getAlignmentFile());
		logToFile("");
		logToFile("MCMC burnIn:\n" + factory.mcmcOptions.burnIn);
		logToFile("MCMC nMCMCSweeps:\n" + factory.mcmcOptions.nMCMCSweeps);
		logToFile("MCMC thinningPeriod:\n" + factory.mcmcOptions.thinningPeriod);
		logToFile("");
		logToFile("Model:\n" + mcmc.model);

		mcmc.run();

		// log running time
		logToFile("Total time in minutes: " + ((System.currentTimeMillis() - startTime)/60000.0));

		// compute the tree
//		MajorityRuleTree.buildAndWriteConsensusTree(
//				Results.getFileInResultFolder("FES.trees.newick"),
//				Results.getFileInResultFolder("FESConsensusTree.Nexus"));
		//
		//		// copy the results to another folder 
		//		File newDirectory = new File(Results.getResultFolder().getParent() + "/experiment." + Results.getResultFolder().getName() + "." + 1 + "." + System.currentTimeMillis());
		//		newDirectory.mkdir();
		//		try {
		//			FileUtils.copyDirectory(Results.getResultFolder(), newDirectory);
		//		} catch (IOException e) {
		//			// TODO Auto-generated catch block
		//			e.printStackTrace();
		//		}
	}

	public MCMCAlgorithm initForwardSampling() {
		factory.addProcessor(this);
		//factory.mcmcOptions.random = new Random(2000);
		factory.mcmcOptions.nMCMCSweeps = nMCMCSweeps;
		factory.mcmcOptions.burnIn = burnIn;
		factory.mcmcOptions.thinningPeriod = thinningPeriod;
		
		// reference bug http://stackoverflow.com/questions/4829576/javac-error-inconvertible-types-with-generics

		//factory.excludeNodeMove((Class<? extends Move>) (Object)SingleBranchScaling.class);
		//factory.excludeNodeMove(AllBranchesScaling.class);
		//factory.excludeNodeMove(PhyloHMCMove.class);
		//factory.excludeNodeMove((Class<? extends Move>) (Object)SingleNNI.class);
		//factory.excludeNodeMove((Class<? extends Move>) (Object)RealVectorMHProposal.class);

		model = new Model();
		MCMCAlgorithm mcmc = factory.build(model, false);

		return mcmc;
	}



	public static void main(String [] args) throws ClassNotFoundException, IOException
	{

		//		MajorityRuleTree.buildAndWriteConsensusTree(
		//				new File
		//				("/Users/sohrab/project/conifer/results/all/2014-07-22-12-03-35-E0EbIo6x.exec/FES.trees.newick"),
		//				new File("/Users/sohrab/project/conifer/results/all/2014-07-22-12-03-35-E0EbIo6x.exec/FESConsensusTree.Nexus"));
		//		return;

		System.out.println("Running the newest version...");

		args = new String[4];
		args[0] = "-initialTreeFilePath";
		args[1] = "/Users/sohrab/project/conifer/src/main/resources/conifer/sampleInput/FES.ape.4.nwk";
		args[2] = "-alignmentFile";
		args[3] = "/Users/sohrab/project/conifer/src/main/resources/conifer/sampleInput/FES_4.fasta";

		args = new String[1];
		args[0] = "-help";

		// TODO: remove this
		for (int i = 0; i < args.length; i++) {
			System.out.println(args[i]);	
		}

		System.out.println(args.length);

		Mains.instrumentedRun(args, new TestPhyloModel());
		
		// TODO: test if all topologies are the same
		// TODO: test if all branch lengths are the same/different 
	}

	public static void makeSyntheticData() {
		List realizations = new ArrayList();
		TestPhyloModel runner = new TestPhyloModel();
		MCMCAlgorithm algo = runner.initForwardSampling();
		ForwardSampler f = new ForwardSampler(algo.model);

		int nSteps = 1;
		Map<Object,List<Double>> results = Maps.newHashMap();
		System.out.println(runner.model.likelihood.observations.toString());
		for (int i = 0; i < nSteps; i++) {
			f.simulate(algo.options.random);
			collectValues(algo, results);
			realizations.add(runner.model.likelihood.observations.toString());
			runner.writeTree(runner.model.likelihood.tree);
			//System.out.println(runner.model.likelihood.observations.toString());
			try {
				FastaUtils.writeFasta(runner.model.likelihood.observations, 
						Results.getFileInResultFolder("SimulatedData.fasta"));
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		TreeObservations observations = runner.model.likelihood.observations;
		// TODO: get the rest of the simulation parameters (total_branch_length?)
		System.out.println("Total BranchLength: " + UnrootedTreeUtils.totalTreeLength((runner.model.likelihood.tree)));
		// TODO: get the tree topology
		System.out.println(results);

	}

	public void writeTree(UnrootedTree tree) {
		PrintWriter theWriter = BriefIO.output(Results.getFileInResultFolder("SimulatedDataTree.newick"));
		theWriter.println(tree.toNewick());
		theWriter.flush();
	}

	public static void collectValues(MCMCAlgorithm mcmc, Map<Object,List<Double>> values) {

		for (Object var: mcmc.model.getLatentVariables()) {
			System.out.println(var.toString());
			System.out.println(mcmc.model.getName(var));
		}

		for (RealValued realValuedVariable : mcmc.model.getLatentVariables(RealValued.class))
			BriefMaps.getOrPutList(values, 
					mcmc.model.getName(realValuedVariable)).add(realValuedVariable.getValue());
		/*
		for (Processor processor : mcmc.processors)
			processor.process(new ProcessorContext(mcmc.options.nMCMCSweeps - 1, mcmc.model, mcmc.options));

		for (RealValued realValuedProcessor : ReflexionUtils.sublistOfGivenType(mcmc.processors, RealValued.class))
			BriefMaps.getOrPutList(values, realValuedProcessor).add(realValuedProcessor.getValue());
		 */
	}

	@Override
	public void process(ProcessorContext context)
	{
		System.out.println("Writing the tree...");
		treeWriter.println(model.likelihood.tree.toNewick());
		treeWriter.flush();
	}

	public void logToFile(String someline) {
		this.detailWriter.println(someline);
		this.detailWriter.flush();
	}

	public String getAlignmentFile() {
		return alignmentFilePath;
	}
}
