package conifer;

/**
 * Created by crystal on 2015-06-07.
 * @author Tingting Zhao (zhaott0416@gmail.com)
 *
 */

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import bayonet.distributions.Normal;
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


public class SimplePhyloSimulatorFixedTreeRateMtx implements Runnable, Processor {
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

    @Option(gloss="File of the tree topology")
    public File treeFile;

    @Option
    public int nTaxa = 5;

    @Option(gloss="select rate matrix model")
    public RateMtxNames selectedRateMtx=RateMtxNames.DNAGTR;


    public List<TreeNode> makeLeaves(int nTaxa, String prefix)
    {
        List<TreeNode> result = new ArrayList<TreeNode>();
        for (int i = 0; i < nTaxa; i++) {
            result.add(TreeNode.withLabel(prefix + i));
        }
        return result;
    }



    public class FixedTopologyAndBranchLengthModel implements SimplePhyloModelContainer
    {
        @DefineFactor(onObservations = true)
        public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood =
                UnrootedTreeLikelihood.createEmptyWithFixedTree(nSites, treeFile, selectedRateMtx)
                        .withExpFamMixture(ExpFamMixture.rateMtxModel(selectedRateMtx));

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

        model = new FixedTopologyAndBranchLengthModel();

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


        try {
            makeSyntheticData(this.treeFile);
        } catch (IOException e) {
            e.printStackTrace();
        }

        // copy the results to another folder
        File newDirectory = new File(
                Results.getResultFolder().getParent() + "/experiment." +
                        Results.getResultFolder().getName() + "." +
                        treeFile + "_"  +
                        System.currentTimeMillis());
        newDirectory.mkdir();
        try {
            FileUtils.copyDirectory(Results.getResultFolder(), newDirectory);
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println(newDirectory.getAbsolutePath());

        //mcmc.run();

        // log running time
        logToFile("Total time in minutes: " + ((System.currentTimeMillis() - startTime)/60000.0));

    }

    public static void main(String [] args)
    {
        Mains.instrumentedRun(args, new SimplePhyloSimulatorFixedTreeRateMtx());
    }

    @SuppressWarnings("unchecked")
    public MCMCAlgorithm getMCMCAlgorithm() {
        factory.addProcessor(this);
        factory.mcmcOptions.nMCMCSweeps = nMCMCSweeps;
        factory.mcmcOptions.burnIn = burnIn;
        factory.mcmcOptions.thinningPeriod = thinningPeriod;
        model = new FixedTopologyAndBranchLengthModel();
        return factory.build(model, false);
    }



    public static void makeSyntheticData(File treeFile) throws IOException {
        //List realizations = new ArrayList();
        SimplePhyloSimulatorFixedTreeRateMtx runner = new SimplePhyloSimulatorFixedTreeRateMtx();
        runner.detailWriter = BriefIO.output(Results.getFileInResultFolder("experiment.details.txt"));

        runner.treeFile = treeFile;

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
            FastaUtils.writeFasta(runner.model.getLikelihood().observations, Results.getFileInResultFolder("SimulatedData.fasta"));
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
