package conifer;

/**
 * Created by crystal on 2015-06-07.
 * @author Tingting Zhao (zhaott0416@gmail.com)
 *
 * The difference between this class and SimplePhyloSimulator is that in this class I can provide the weights
 * for both univariate features and bivariate features in a json file in the format of a json map and we can
 * generate the sequences under a provided fixed tree and fixed weights values corresponding to specific features
 * in the json file so that we can generate data according to a specific rate matrix. We need to provide the json
 * file as a command line option for weightFile
 */

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.HashMap;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;

import bayonet.distributions.Normal;
import conifer.ctmc.expfam.ExpFamParameters;
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
    public RateMtxNames selectedRateMtx=RateMtxNames.POLARITYSIZEGTR;

    @Option(gloss="Weights values for the rate matrix")
    public File weightFile= new File("/Users/crystal/Documents/codeExperiment/Revision/Experiment2/RealData/RealWeights.txt");

    public FileInputStream in=null;

    //public double [] weights;
    public Map<String, Double> featureWeights = new HashMap<String, Double>();

    public List<TreeNode> makeLeaves(int nTaxa, String prefix)
    {
        List<TreeNode> result = new ArrayList<TreeNode>();
        for (int i = 0; i < nTaxa; i++) {
            result.add(TreeNode.withLabel(prefix + i));
        }
        return result;
    }

    public FileInputStream setFileInputStream(File file){
        try {
            in = new FileInputStream(file);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return in;
    }

    /**
     * This method will read the provided json map  "file" to create a HashMap with keys as the features
     * to create the rate matrix and the values for the map are the weights values corresponding to the
     * features.
     * @param file
     * @return
     */

    public Map<String, Double> setMapFeatureWeights(File file){

        ObjectMapper mapper = new ObjectMapper();
        FileInputStream curIn = setFileInputStream(file);

        try {

            TypeReference<HashMap<String,Double>> typeRef
                    = new TypeReference<HashMap<String,Double>>(){};

            featureWeights = mapper.readValue(curIn, typeRef);

            //weights = mapper.readValue(curIn, double [].class);
        } catch (IOException e) {
            e.printStackTrace();
        }finally{
            if (curIn != null)
                try {
                    curIn.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
        }

        return featureWeights;
    }

    public class FixedTopologyAndBranchLengthModel implements SimplePhyloModelContainer
    {

        @DefineFactor(onObservations = true)
        public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood =
                UnrootedTreeLikelihood.createEmptyWithFixedTree(nSites, treeFile, selectedRateMtx)
                        .withExpFamMixture(ExpFamMixture.rateMtxModel(selectedRateMtx));



        @DefineFactor(onObservations = true)
        public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior =
                IIDRealVectorGenerativeFactor
                        .iidNormalOn(likelihood.evolutionaryModel.rateMatrixMixture.parameters);

        public UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> getLikelihood() {
            return this.likelihood;
        }

        /**
         * this method will help us set the weight values when creating the rate matrix. It will set
         * the values for specific features based on the hashmap featureWeights
         * @param featureWeights
         */
        public void setVectorWeights(Map<String, Double> featureWeights) {
            likelihood.evolutionaryModel.rateMatrixMixture.parameters.setVector(featureWeights);
        }
    }

    // we can build different models here and instantiate them according to the input parameter
    // in the run block.

    // Note: only instantiate this in run to avoid problems with command line argument parsing
    public FixedTopologyAndBranchLengthModel model;

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

        setMapFeatureWeights(weightFile);
        model = new FixedTopologyAndBranchLengthModel();
        model.setVectorWeights(featureWeights);

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



    public void makeSyntheticData(File treeFile) throws IOException {
        //List realizations = new ArrayList();
        SimplePhyloSimulatorFixedTreeRateMtx runner = new SimplePhyloSimulatorFixedTreeRateMtx();
        runner.detailWriter = BriefIO.output(Results.getFileInResultFolder("experiment.details.txt"));

        runner.treeFile = treeFile;

        runner.featureWeights= runner.setMapFeatureWeights(weightFile);

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
