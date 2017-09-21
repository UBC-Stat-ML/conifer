package conifer;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;

import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.ForwardSampler;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.mcmc.Move;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import blang.variables.RealValued;
import briefj.BriefIO;
import briefj.BriefMaps;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;

import com.google.common.collect.Maps;

import conifer.TestPhyloModel.Model;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.RateMtxNames;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.io.FastaUtils;
import conifer.io.TreeObservations;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.moves.SPRMove;
import conifer.moves.SingleNNI;


/**
 * This class will simulate data from TestPhyloModel
 * @author Sohrab Salehi (sohrab.salehi@gmail.com)
 *
 */
public class SimplePhyloSimulator implements Runnable, Processor {
    @Option
	public int nMCMCSweeps = 10000;

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

	@Option
	public int nTaxa = 5;
	
	@Option
	public boolean fixedTopology = false;
	
	@Option
	public boolean fixedBranchLength = false;
	
	@Option(gloss="select rate matrix model")
	public static RateMtxNames selectedRateMtx=RateMtxNames.DNAGTR;

	@Option(gloss="Indicator of we normalize the rate matrix if it is set to true")
	public boolean isNormalized = true;
	
	public List<TreeNode> makeLeaves(int nTaxa, String prefix) 
	{	
		List<TreeNode> result = new ArrayList<TreeNode>();
		for (int i = 0; i < nTaxa; i++) {
			result.add(TreeNode.withLabel(prefix + i));
		}
		return result;
	}
	
	public class FullLatentModel implements SimplePhyloModelContainer
	{
		@DefineFactor(onObservations = true)
		public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = 
		UnrootedTreeLikelihood.createEmpty(nSites, makeLeaves(nTaxa, "t"), selectedRateMtx)		
		.withExpFamMixture(ExpFamMixture.rateMtxModel(selectedRateMtx, isNormalized));

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
		
		public UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> getLikelihood() {
			return this.likelihood;
		}
	}

	
	public class FixedTopologyAndBranchLengthModel implements SimplePhyloModelContainer
	{
		@DefineFactor(onObservations = true)
		public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = 
		UnrootedTreeLikelihood.createEmpty(nSites, makeLeaves(nSites, "t"), selectedRateMtx)		
		.withExpFamMixture(ExpFamMixture.rateMtxModel(selectedRateMtx, isNormalized));

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
		
		public UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> getLikelihood() {
			return this.likelihood;
		}
	}
	
	// we can build different models here and instantiate them according to the input parameter
	// in the run block.

	// Note: only instantiate this in run to avoid problems with command line argument parsing
	public SimplePhyloModelContainer model;

	private PrintWriter treeWriter;
	private PrintWriter detailWriter;

	@SuppressWarnings("unchecked")
	@Override
	public void run()
	{
		System.out.println("Running the newest version...");
		treeWriter = BriefIO.output(Results.getFileInResultFolder("FES.trees.newick"));
		detailWriter = BriefIO.output(Results.getFileInResultFolder("experiment.details.txt"));
		
		factory.addProcessor(this);
		factory.mcmcOptions.nMCMCSweeps = nMCMCSweeps;
		factory.mcmcOptions.burnIn = burnIn;
		factory.mcmcOptions.CODA = true;
		factory.mcmcOptions.thinningPeriod = thinningPeriod;
		
		if (fixedTopology && fixedBranchLength) {
			model = new FixedTopologyAndBranchLengthModel();
		} else if (fixedTopology && !fixedBranchLength) {
			model = new FullLatentModel();
			factory.excludeNodeMove((Class<? extends Move>) (Object)SingleNNI.class);
			factory.excludeNodeMove(SPRMove.class);
		} else {
			model = new FullLatentModel();
		}
		
		MCMCAlgorithm mcmc = factory.build(model, false);
		System.out.println(mcmc.model);

		long startTime = System.currentTimeMillis();
		String excluding = "GTR-Excluding all but SPR move";

		// log experiment information
		logToFile("Experiment Title:" + excluding);
		logToFile("Current Date:" + RunFacility.getCurrentDateString());
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
		MajorityRuleTree.buildAndWriteConsensusTree(
				Results.getFileInResultFolder("FES.trees.newick"),
				Results.getFileInResultFolder("FESConsensusTree.Nexus"));
	}

	@SuppressWarnings("unchecked")
	public MCMCAlgorithm getMCMCAlgorithm() {
		factory.addProcessor(this);
		factory.mcmcOptions.nMCMCSweeps = nMCMCSweeps;
		factory.mcmcOptions.burnIn = burnIn;
		factory.mcmcOptions.thinningPeriod = thinningPeriod;

		if (fixedTopology && fixedBranchLength) {
			model = new FixedTopologyAndBranchLengthModel();
		} else if (fixedTopology && !fixedBranchLength) {
			model = new FullLatentModel();
			factory.excludeNodeMove((Class<? extends Move>) (Object)SingleNNI.class);
			factory.excludeNodeMove(SPRMove.class);
		} else {
			model = new FullLatentModel();
		}
		
		return factory.build(model, false);
	}


	public static void main(String [] args) throws ClassNotFoundException, IOException
	{
		int[] numbersOfTaxa = {5, 10, 50};

		for (int i : numbersOfTaxa) {
			for (ComparisonModes mode :  ComparisonModes.values()) {
				makeSyntheticData(i, mode, selectedRateMtx);
				
				// copy the results to another folder 
				File newDirectory = new File(
						Results.getResultFolder().getParent() + "/experiment." + 
				Results.getResultFolder().getName() + "." + 
								i + "_" + mode.toString() + "_" +
				System.currentTimeMillis());
				newDirectory.mkdir();
				FileUtils.copyDirectory(Results.getResultFolder(), newDirectory);
				System.out.println(newDirectory.getAbsolutePath());
			}
		}
	}

	public static void makeSyntheticData(int numberOfTaxa, ComparisonModes mode, RateMtxNames selectedRateMtx) throws IOException {
		//List realizations = new ArrayList();
		SimplePhyloSimulator runner = new SimplePhyloSimulator();
		runner.detailWriter = BriefIO.output(Results.getFileInResultFolder("experiment.details.txt"));
		
		runner.nTaxa = numberOfTaxa;
		
		switch (mode) {
		case DEFAULT:
			runner.fixedBranchLength = false;
			runner.fixedTopology = false;
			break;
		case FIXED_TOPOLOGY_LATENT_BRANCH_LENGTHS:
			runner.fixedBranchLength = false;
			runner.fixedTopology = true;
			break;
		case FIXED_TOPOLOGY_FIXED_BRANCH_LENGHTS:
			runner.fixedBranchLength = true;
			runner.fixedTopology = true;
			break;
		default:
			System.out.println("Illigal mode!");
			break;
		}
		
		MCMCAlgorithm algo = runner.getMCMCAlgorithm();
		ForwardSampler f = new ForwardSampler(algo.model);

		int nSteps = 1;
		Map<Object,List<Double>> results = Maps.newHashMap();
		System.out.println(runner.model.getLikelihood().observations.toString());
		for (int i = 0; i < nSteps; i++) {
			f.simulate(algo.options.random);
			//collectValues(algo, results);
			//realizations.add(runner.model.getLikelihood().observations.toString());
			runner.writeTree(runner.model.getLikelihood().tree);
			// write the FASTA file corresponding to the simulation
			FastaUtils.writeFasta(runner.model.getLikelihood().observations, Results.getFileInResultFolder("SimulatedData.fasta"), selectedRateMtx);
		}
		
		runner.detailWriter.write("Total BranchLength: " +
		UnrootedTreeUtils.totalTreeLength((runner.model.getLikelihood().tree)) + "\n");

		// TODO: get the rest of the simulation parameters (total_branch_length?)
		for (Object var: algo.model.getLatentVariables()) {
			runner.detailWriter.write(algo.model.getName(var) + " : " + var.toString() + "\n");
			System.out.println(algo.model.getName(var) + " : " + var.toString());
		}
		
		runner.detailWriter.flush();
	}

	public void writeTree(UnrootedTree tree) {
		PrintWriter theWriter = BriefIO.output(Results.getFileInResultFolder("SimulatedDataTree.newick"));
		theWriter.println(tree.toNewick());
		theWriter.flush();
	}

	public static void collectValues(MCMCAlgorithm mcmc, Map<Object,List<Double>> values) {

//		for (Object var: mcmc.model.getLatentVariables()) {
//			System.out.println(var.toString());
//			System.out.println(mcmc.model.getName(var));
//		}

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
		treeWriter.println(model.getLikelihood().tree.toNewick());
		treeWriter.flush();
	}

	public void logToFile(String someline) {
		this.detailWriter.println(someline);
		this.detailWriter.flush();
	}


}
