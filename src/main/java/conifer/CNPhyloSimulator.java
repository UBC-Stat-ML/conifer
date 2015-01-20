package conifer;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;

import blang.ForwardSampler;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.annotations.DefineFactor;
import blang.mcmc.Move;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import briefj.BriefIO;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;
import conifer.ctmc.cnv.CopyNumberMixture;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.io.CopyNumberTreeObservation;
import conifer.models.CNMultiCategorySubstitutionModel;
import conifer.models.CNPair;
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
	@Option
	public final int nSites = 100;
    @Option
	public final int nTaxa = 3;
	
	
	public class Model {
	
		
		Set<TreeNode> leaves = new HashSet<TreeNode>(TopologyUtils.makeLeaves(nTaxa, "t_"));
		
		@DefineFactor
		public final UnrootedTreeLikelihood<CNMultiCategorySubstitutionModel<CopyNumberMixture>> likelihood = UnrootedTreeLikelihood.createEmptyCN(nSites, leaves);

	      
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
		
			// FastaUtils.writeFasta(model.likelihood.observations, Results.getFileInResultFolder("SimulatedData.fasta"));
			logObservations( (CopyNumberTreeObservation) model.likelihood.observations);
			// TODO: get the rest of the simulation parameters, if any
			for (Object var: mcmc.model.getLatentVariables()) {
				System.out.println(mcmc.model.getName(var) + " : " + var.toString());
			}
		}
	}
	
	
	public void logObservations(CopyNumberTreeObservation observation) {
		System.out.println(observation.toString());
		
		// Write CTMC State, that is the copy-numbers
		writeCTMC(observation);
		
	}
	
	public void writeCTMC(CopyNumberTreeObservation observation) {
		StringBuilder result = new StringBuilder();
		
		// Add the header
		result.append("\"sample_id\",\"site_id\",\"ref_counts\",\"alt_counts\",\"cluster_id\",\"ref_CN\",\"alt_CN\"" + "\n");
		
		Set<TreeNode> leaves = observation.getLeaves();
		
		LinkedHashMap<TreeNode, List<CNPair>> cns = observation.getPrintFriendlyCTMCState();
		LinkedHashMap<TreeNode, List<CNPair>> emissions = observation.getEmissions();
		
	    for (TreeNode node : leaves) {
	        // do not record root level
	        if (node.toString() != "root")
	        {
	        List<CNPair> cnPairs = cns.get(node);
	    	for (int i = 0; i < cnPairs.size(); i++) {
	    		CNPair currPair = cnPairs.get(i);
	    		result.append(node.toString() +  ", " + // add sample_id
	    		    		  i + ", " + 				// site_id
	    		    		  emissions.get(node).get(i).getrA() + "," + // emission rA
	    		    		  emissions.get(node).get(i).getRa() + "," + // emission ra
	    		    		  node.toString() + ","  +                // cluster == one per sample
	    		    		  currPair.getrA() + ", " + // ref_CN
	    		    		  currPair.getRa() + ", " + // alt_CN
	    		    		  "\n");
			}
	        }
	    }
	    
	    // write the result to file
	    PrintWriter theWriter = BriefIO.output(Results.getFileInResultFolder("simulatedCNData.csv"));
		theWriter.println(result.toString());
		theWriter.flush();
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
