package conifer;

import java.io.File;
import java.io.IOException;
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
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;


public class CNPhyloModel extends MCMCRunner
{

    @Option(required = true, gloss = "file containing raw reads") String emissionData;
    
    File inputFile = new File(emissionData);

    @DefineFactor(onObservations = true)
    public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = 
    UnrootedTreeLikelihood
    .fromFastaFile(inputFile)
    .withExpFamMixture(ExpFamMixture.kimura1980());
   

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
        CNPhyloModel runner = new CNPhyloModel();
        runner.factory.mcmcOptions.nMCMCSweeps = 10000;
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
