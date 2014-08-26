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
import briefj.run.Results;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;



public class SimplePhyloModel extends MCMCRunner
{
	File inputFile 
	= new File("src/main/resources/conifer/sampleInput/FES_4.fasta");

	@DefineFactor(onObservations = true)
	public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = 
	UnrootedTreeLikelihood
	.fromFastaFile(inputFile)
	.withExpFamMixture(ExpFamMixture.kimura1980())
	.withTree(new File("src/main/resources/conifer/sampleInput/FES.ape.4.nwk"));

	@DefineFactor
	NonClockTreePrior<RateParameterization> treePrior = 
	NonClockTreePrior
	.on(likelihood.tree);

	@DefineFactor
	Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = 
	Exponential
	.on(treePrior.branchDistributionParameters.rate)
	.withMean(10.0);

	@DefineFactor
	public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior =
	IIDRealVectorGenerativeFactor
	.iidNormalOn(likelihood.evolutionaryModel.rateMatrixMixture.parameters);

	private final PrintWriter treeWriter = BriefIO.output(Results.getFileInResultFolder("FES.trees.newick"));

	public static void main(String [] args) throws ClassNotFoundException, IOException
	{
		SimplePhyloModel runner = new SimplePhyloModel();
		runner.factory.mcmcOptions.nMCMCSweeps = 100;
		runner.factory.mcmcOptions.burnIn = (int) Math.round(.1 * runner.factory.mcmcOptions.nMCMCSweeps);
		
		// run
		runner.run();

		// compute the tree
		MajorityRuleTree.buildAndWriteConsensusTree(
				Results.getFileInResultFolder("FES.trees.newick"),
				Results.getFileInResultFolder("FESConsensusTree.Nexus"));
	}

	protected void process(ProcessorContext context)
	{
		treeWriter.println(likelihood.tree.toNewick());
		treeWriter.flush();
	}

	public String getInputFile() {
		return inputFile.getAbsolutePath();
	}
}
