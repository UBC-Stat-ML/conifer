package conifer;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.RateMtxNames;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.mcmc.Move;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import briefj.BriefIO;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.moves.AllBranchesScaling;
import conifer.moves.PhyloHMCMove;
import conifer.moves.RealVectorMHProposal;
import conifer.moves.SPRMove;
import conifer.moves.SingleBranchScaling;
import conifer.moves.SingleNNI;

public class InstrumentedSimplePhyloModel implements Runnable, Processor {
	@Option
	public String initialTreeFilePath = 
	new String("src/main/resources/conifer/sampleInput/FES.ape.4.nwk");

	@Option
	public String alignmentFilePath =  
	new String("src/main/resources/conifer/sampleInput/FES_4.fasta");

	@Option
	public int nMCMCSweeps = 100;

	@Option
	public int burnIn = (int) Math.round(.1 * nMCMCSweeps);

	@Option
	public int thinningPeriod = 10;

	@OptionSet(name = "factory")
	public final MCMCFactory factory = new MCMCFactory();


	@Option
	public int nSites = 500;
	
	@Option(gloss="provide rate matrix method") 
	public RateMtxNames selectedRateMtx;

	@OptionSet(name = "NodeMoves")
	public final NodeMoves nd = new NodeMoves();

	@Option(gloss="Indicator of we normalize the rate matrix if it is set to true")
	public boolean isNormalized = true;
	
    public static class NodeMoves
    {
        @Option public boolean SPRMove = true; 
        @Option public boolean SingleNNI = true;
        @Option public boolean SingleBranchScaling = true;
        @Option public boolean AllBranchesScaling = true;
        @Option public boolean PhyloHMCMove = true;
        @Option public boolean RealVectorMHProposal = true;         
    }
	
	 
	
	public class Model
	{
		@DefineFactor(onObservations = true)
		public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = 
		UnrootedTreeLikelihood
		.fromFastaFile(new File(alignmentFilePath), selectedRateMtx.KIMURA1980)
		.withExpFamMixture(ExpFamMixture.rateMtxModel(selectedRateMtx,isNormalized))
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


	// Note: only instantiate this in run to avoid problems with command line argument parsing
	public Model model;

	private PrintWriter treeWriter;

	@Override
	public void run()
	{
		treeWriter = BriefIO.output(Results.getFileInResultFolder("FES.trees.newick"));

		factory.addProcessor(this);
		factory.mcmcOptions.nMCMCSweeps = nMCMCSweeps;
		factory.mcmcOptions.burnIn = burnIn;
		factory.mcmcOptions.CODA = true;
		factory.mcmcOptions.thinningPeriod = thinningPeriod;

		if(!nd.SPRMove)
		    factory.excludeNodeMove(SPRMove.class); 
        if(!nd.AllBranchesScaling)
            factory.excludeNodeMove(AllBranchesScaling.class);
        if(!nd.PhyloHMCMove)
            factory.excludeNodeMove(PhyloHMCMove.class);        
        if(!nd.RealVectorMHProposal)
            factory.excludeNodeMove(RealVectorMHProposal.class); 		
        if(!nd.SingleBranchScaling)
            factory.excludeNodeMove(SingleBranchScaling.class);
        if(!nd.SingleNNI)
            factory.excludeNodeMove(SingleNNI.class);
		
		// init the model
		model = new Model();
		MCMCAlgorithm mcmc = factory.build(model, false);
		
        File graph = Results.getFileInResultFolder("probability-graph.dot");  
        mcmc.model.printGraph(graph);
		
		mcmc.run();

		// compute the tree
		MajorityRuleTree.buildAndWriteConsensusTree(
				Results.getFileInResultFolder("FES.trees.newick"),
				Results.getFileInResultFolder("FESConsensusTree.Nexus"));
	}


	public static void main(String [] args) throws ClassNotFoundException, IOException
	{		
		Mains.instrumentedRun(args, new InstrumentedSimplePhyloModel());
	}

	
	public void writeTree(UnrootedTree tree) {
		PrintWriter theWriter = BriefIO.output(Results.getFileInResultFolder("SimulatedDataTree.newick"));
		theWriter.println(tree.toNewick());
		theWriter.flush();
	}

	@Override
	public void process(ProcessorContext context)
	{
		treeWriter.println(model.likelihood.tree.toNewick());
		treeWriter.flush();
	}

	public String getAlignmentFile() {
		return alignmentFilePath;
	}

}
