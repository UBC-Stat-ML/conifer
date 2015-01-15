package conifer;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;

import bayonet.distributions.Uniform;
import bayonet.distributions.Uniform.MinMaxParameterization;
import blang.ForwardSampler;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.annotations.DefineFactor;
import blang.mcmc.Move;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import briefj.BriefIO;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;
import conifer.ctmc.cnv.CopyNumberMixture;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.moves.SPRMove;
import conifer.moves.SingleNNI;


/*
 * TODO: fix the tree, only simulate copy number change and one eventual point-mutation,
 * TODO: later, add the option to jump between trees
 * 
 */

public class CNPhyloSimulator implements Runnable, Processor {
	

	@OptionSet(name = "factory")
	public final MCMCFactory factory = new MCMCFactory();

	
	public class Model {

		int nSites = 5;
		int nTaxa = 5;
		
		Set<TreeNode> leaves = new HashSet<TreeNode>(TopologyUtils.makeLeaves(nTaxa, "t_"));
		
		@DefineFactor
		public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<CopyNumberMixture>> likelihood = UnrootedTreeLikelihood.createEmptyCN(nSites, leaves);
		
	     @DefineFactor
	     Uniform<MinMaxParameterization> cnvParameter = Uniform.on(likelihood.evolutionaryModel.rateMatrixMixture.parameters.alpha).withBounds(0.25, 0.75); 
	      
	}

	public Model model;
	
	@SuppressWarnings("unchecked")
	@Override
	public void run() 
	{
		model = new Model();
		MCMCAlgorithm mcmc = factory.build(model, false);
		
		// Fix topology
		factory.excludeNodeMove((Class<? extends Move>) (Object) SingleNNI.class);
		factory.excludeNodeMove(SPRMove.class);
		
		ForwardSampler f = new ForwardSampler(mcmc.model);
		
		// Simulate K data sets
		int nDataSets = 1;
		for (int i = 0; i < nDataSets; i++) {
			f.simulate(mcmc.options.random);
			
			// write the CN-input file corresponding to the simulation
			System.out.println(model.likelihood.observations.toString());
			// FastaUtils.writeFasta(model.likelihood.observations, Results.getFileInResultFolder("SimulatedData.fasta"));

			// TODO: get the rest of the simulation parameters, if any
			for (Object var: mcmc.model.getLatentVariables()) {
				System.out.println(mcmc.model.getName(var) + " : " + var.toString());
			}
		}
	}
	
	public void writeTree(UnrootedTree tree) {
		PrintWriter theWriter = BriefIO.output(Results.getFileInResultFolder("SimulatedDataTree.newick"));
		theWriter.println(tree.toNewick());
		theWriter.flush();
	}
	
	public static void main(String[] args) throws ClassNotFoundException, IOException {
		Mains.instrumentedRun(args, new CNPhyloSimulator());
	}


	@Override
	public void process(ProcessorContext context) {
		// TODO Auto-generated method stub
		
	}
}
