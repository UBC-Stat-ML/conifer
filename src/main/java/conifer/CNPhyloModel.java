package conifer;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.Uniform;
import bayonet.distributions.Uniform.MinMaxParameterization;
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
import conifer.models.MultiCategorySubstitutionModel;
import conifer.moves.AllBranchesScaling;
import conifer.moves.SPRMove;

public class CNPhyloModel implements Runnable, Processor 
{
    @Option(required = true, gloss = "file containing raw reads")
    public String emissionData;

    @OptionSet(name = "factory")
    public final MCMCFactory factory = new MCMCFactory();

    public class Model 
    {
        File inputFile = new File(emissionData);

        @DefineFactor
        public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<CopyNumberMixture>> likelihood = UnrootedTreeLikelihood
        .fromCNFile(inputFile);

        @DefineFactor
        NonClockTreePrior<RateParameterization> treePrior = NonClockTreePrior.on(likelihood.tree);

        @DefineFactor
        Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = Exponential.on(
                treePrior.branchDistributionParameters.rate).withMean(10.0);

//        @DefineFactor
//        Uniform<MinMaxParameterization> cnvParameter = Uniform.on(likelihood.evolutionaryModel.rateMatrixMixture.parameters.alpha).withBounds(0.25, 0.75); 

    }

    public Model model;

    private PrintWriter treeWriter;

    @Override
    public void run() 
    {
        treeWriter = BriefIO.output(Results.getFileInResultFolder("cntrees.newick"));
        factory.addProcessor(this);
        factory.excludeNodeMove(SPRMove.class);
        factory.excludeNodeMove(AllBranchesScaling.class);
        model = new Model();
        MCMCAlgorithm mcmc = factory.build(model, false);
        File graph = Results.getFileInResultFolder("probability-graph.dot");	
        mcmc.model.printGraph(graph);
        mcmc.options.CODA = true;
        System.out.println(mcmc.model.toString());

        mcmc.run();

        MajorityRuleTree.buildAndWriteConsensusTree(Results.getFileInResultFolder("cntrees.newick"),
                Results.getFileInResultFolder("cnConsensusTree.Nexus"));

    }

    public static void main(String[] args) throws ClassNotFoundException, IOException
    {
        Mains.instrumentedRun(args, new CNPhyloModel());
    }

    @Override
    public void process(ProcessorContext context)
    {
        treeWriter.println(model.likelihood.tree.toNewick());
        treeWriter.flush();
    }

}
