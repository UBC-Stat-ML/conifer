package conifer;

import java.io.File;
import java.io.PrintWriter;
import java.util.Random;

import conifer.ctmc.EigenCTMC;
import conifer.ctmc.EndPointSampler;
import conifer.moves.RealVectorOverRelaxedSlice;
import org.apache.commons.io.FileUtils;

import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.ForwardSampler;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.MCMCRunner;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.processing.Processor;
import briefj.opt.Option;
import briefj.opt.OptionsParser;
import briefj.opt.OptionSet;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.RateMtxNames;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.moves.PhyloHMCMove;
import conifer.moves.RealVectorMHProposal;
import conifer.processors.FileNameString;
import briefj.run.Mains;
import briefj.run.Results;
import blang.processing.ProcessorContext;
import briefj.BriefIO;
import briefj.OutputManager;

import java.io.IOException;



public class SingleProteinModel implements Runnable, Processor
{
    @Option(gloss="File of provided alignment")
    public File inputFile;


    @Option(gloss="File of the tree topology")
    public File treeFile;

    @Option(gloss="Indicator of whether to exclude HMC move")
    public boolean isExcluded=false;

    @Option(gloss="Indicator of use slice sampler or not")
    public boolean sliceSampler=false;

    @Option(gloss="Number of MCMC iterations")
    public int nMCMCIterations = 100000;

    @Option(gloss="ESS Experiment Number")
    public int rep = 2;

    @Option(gloss="Rate Matrix Method")
    public RateMtxNames selectedRateMtx = RateMtxNames.PAIR;

    @OptionSet(name = "factory")
    public final MCMCFactory factory = new MCMCFactory();

    @Option(gloss="Bandwidth of proposal for vectors in MCMC")
    public double bandwidth = 0.01;

    @Option(gloss="Epsilon provided ")
    public static Double epsilon = null;

    @Option(gloss="L provided")
    public static Integer L = null;

    @Option(gloss="provided size of adaptation")
    public static Integer sizeAdapt = 500;

    @Option(gloss="provided number of HMC iterations per refreshment of Z")
    public static int nItersPerPathAuxVar = 100;

    @Option(gloss="Use cache or not")
    public boolean cached=true;

    @Option(gloss="record cache size or not")
    public boolean recordCacheSize=true;

    @Option(gloss ="Indicator of using Diagonalization of matrix exponential algorithm or not")
    public boolean useDiag = false;

    @Option(gloss="Indicator of using auxiliary variable Z or not when getting the expected statistics of holding time and number of transitions")
    public boolean useAuxiliaryVariable = true;

    public class Model
    {
        @DefineFactor(onObservations = true)
        public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood1 =
                UnrootedTreeLikelihood
                        .fromFastaFile(inputFile, selectedRateMtx)
                        .withExpFamMixture(ExpFamMixture.rateMtxModel(selectedRateMtx))
                        .withTree(treeFile);

        @DefineFactor
        public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior =
                IIDRealVectorGenerativeFactor
                        .iidNormalOn(likelihood1.evolutionaryModel.rateMatrixMixture.parameters);
    }

    private final PrintWriter detailWriter = BriefIO.output(Results.getFileInResultFolder("experiment.details.txt"));

// The topology of the tree is fixed so that I don't put a prior on the tree topology
// @DefineFactor
// NonClockTreePrior<RateParameterization> treePrior1 = 
//  NonClockTreePrior
// .on(likelihood1.tree);


//  @DefineFactor
//  Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = 
//    Exponential
//    .on(treePrior.branchDistributionParameters.rate)
//    .withMean(10.0);

// private final PrintWriter treeWriter = BriefIO.output(Results.getFileInResultFolder("tree.nwk"));

    // Note: only instantiate this in run to avoid problems with command line argument parsing
    public Model model;

    public void run()
    {
        RealVectorMHProposal.bandWidth = bandwidth;
        PhyloHMCMove.epsilon=epsilon;
        PhyloHMCMove.L = L;
        PhyloHMCMove.sizeAdapt=sizeAdapt;
        PhyloHMCMove.nItersPerPathAuxVar=nItersPerPathAuxVar;
        PhyloHMCMove.useAuxiliaryVariable = useAuxiliaryVariable;
        EndPointSampler.cached=cached;
        MultiCategorySubstitutionModel.recordCacheSize=recordCacheSize;
        EigenCTMC.useDiag = useDiag;
        factory.addProcessor(this);
        model = new Model();
        long startTime = System.currentTimeMillis();
        if(isExcluded)
        {
            if(sliceSampler){
                throw new RuntimeException();
            }else{
                // the default setting for isExcluded should be false
                factory.excludeNodeMove(PhyloHMCMove.class);
            }

        }else{

            if(sliceSampler){}else{
                factory.excludeNodeMove(RealVectorOverRelaxedSlice.class);

            }

        }

        if(sliceSampler){
            // when sliceSampler is true, we need to make sure isExcluded is false
            if(isExcluded)
                throw new RuntimeException();
            else{
                factory.excludeNodeMove(PhyloHMCMove.class);
                factory.excludeNodeMove(RealVectorMHProposal.class);

            }

        }

        MCMCAlgorithm mcmc = factory.build(model, false);
        mcmc.options.random = new Random(rep);
        mcmc.options.nMCMCSweeps = nMCMCIterations;
        mcmc.options.thinningPeriod=factory.mcmcOptions.thinningPeriod;
        mcmc.options.burnIn = (int) Math.round(.1 * factory.mcmcOptions.nMCMCSweeps);
        mcmc.run();
        String fileName = inputFile.getName();
        FileNameString fileNameString = new FileNameString(fileName);
        String numberOfSites = fileNameString.subStringBetween(fileName, "numSites", "Seed");
        String whichSeedUsed = fileNameString.subStringBetween(fileName, "Seed", ".txt");
        logToFile("Total time in minutes: " + ((System.currentTimeMillis() - startTime)/60000.0));
        //File newDirectory = new File(Results.getResultFolder().getParent() + "rep"+ rep+ "isExcludedHMCMove" + isExcluded + bandwidth+selectedRateMtx+"numSites"+numberOfSites+"Seed"+whichSeedUsed+ "epsilon"+PhyloHMCMove.epsilon+"L"+PhyloHMCMove.L);
        File newDirectory = new File(Results.getResultFolder().getParent() + "rep"+ rep+ "isExcludedHMCMove" + isExcluded + "slice"+sliceSampler+bandwidth+selectedRateMtx+"numSites"+numberOfSites+"Seed"+whichSeedUsed +"epsilon"+PhyloHMCMove.epsilon+"L"+PhyloHMCMove.L+"Adapt"+ PhyloHMCMove.sizeAdapt+
                "nItersPerPathAuxVar"+nItersPerPathAuxVar);
        newDirectory.mkdir();
        try
        {
            FileUtils.copyDirectory(Results.getResultFolder(), newDirectory);
        } catch (IOException e)
        {
            e.printStackTrace();
        }

    }

    public static void main(String [] args)
    {
        Mains.instrumentedRun(args, new SingleProteinModel());
    }

    @Override
    public void process(ProcessorContext context)
    {

    }
    public void logToFile(String someline) {
        this.detailWriter.println(someline);
        this.detailWriter.flush();
    }


}
