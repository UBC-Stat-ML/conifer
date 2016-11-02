package conifer;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Map;

import briefj.BriefIO;
import briefj.run.Results;
import com.beust.jcommander.internal.Maps;

import conifer.ctmc.expfam.RateMtxNames;
import conifer.global.MultipleConvexCollisionSolver;
import conifer.global.SelfImplementedSolver;
import conifer.global.SelfImplementedSolverNames;
import conifer.global.SolverNames;
import conifer.moves.RealVectorAdaptiveMHProposal;
import conifer.moves.RealVectorOverRelaxedSlice;
import conifer.rejfreemodels.phylo.PhyloLocalRFMove;
import conifer.rejfreemodels.phylo.PhyloRFMove;
import conifer.rejfreeutil.RFSamplerOptions;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.mcmc.Move;
import blang.processing.LogDensityProcessor;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import briefj.opt.InputFile;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.ExpFamParameters;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.moves.PhyloHMCMove;
import conifer.moves.RealVectorMHProposal;
import org.apache.commons.io.FileUtils;


/**
 * Test the phylogenetic MCMC moves on a phylogenetic model.
 *
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 */
public class TestSimuDataLocalSampler implements Runnable, Processor {
    @InputFile
    @Option(required = true, gloss = "Location of sequence files in FASTA format.")
    public File treeFile;

    @InputFile
    @Option(required = true, gloss = "Location of tree file in newick format")
    public File sequencesFile;

    @Option(gloss = "If the Hamiltonian Monte Carlo sampler should be used.")
    public boolean useHMC = false;

    @Option(gloss = "If the random walk (symmetric proposal) Metropolis sampler should be used.")
    public boolean useMH = false;

    @Option(gloss = "If L and epsilon should be set adaptively.")
    public boolean useAdaptiveHMC = false;

    @Option(gloss = "Step size of the the HMC sampler")
    public double epsilon = 0.05;

    @Option(gloss = "Number of steps per accept reject HMC move.")
    public int L = 100;

    @Option(gloss = "Number of parameter moves per auxiliary variable resampling.")
    public int nItersPerPathAuxVar = 1000;

    @Option(gloss = "If the Rejection Free sampler should be used.")
    public boolean useGlobalRF = true;

    @Option(gloss = "If the local Rejection Free Sampler should be used")
    public boolean useLocalRF = false;

    @Option(gloss="Indicator of we normalize the rate matrix if it is set to true")
    public boolean isNormalized = false;

    @Option(gloss = "If the local rejection free sampler is used, we provide the fixed trajectory length")
    public Double maxTrajectoryLength=null;


    @Option(gloss="Number of MCMC iterations")
    public int nMCMCIterations = 1;


    @Option(gloss="Rate Matrix Method")
    public RateMtxNames selectedRateMtx = RateMtxNames.PROTEINSIMPLEGTR;

    @Option(gloss="Use superposition method in global sampler")
    public boolean useSuperPosition = false;

    @Option(gloss="Use Pegasus solver in global sampler or not, if not use Brent Solver")
    public boolean usePegasusSolver = false;

    @Option(gloss="Use self implemented solver or not, if not use Pegasus Solver")
    public boolean useSelfImplementedSolver = true;

    @Option(gloss = "self implemented solver names")
    public SelfImplementedSolverNames selfSolverNames = SelfImplementedSolverNames.ImprovedPegasus;

    @OptionSet(name = "rfoptions")
    public RFSamplerOptions rfOptions = new RFSamplerOptions();

    @OptionSet(name = "factory")
    public final MCMCFactory factory = new MCMCFactory();

    @Option(gloss = "lowerbound of the trajectory length for local rejection free sampler")
    public double lowerbound = 0.02;

    @Option(gloss =  "upperbound of the trajectory length for local rejection free sampler")
    public double upperbound = 1.0;

    @Option(gloss = "If the Adaptive Local Rejection Free sampler should be used.")
    public boolean useAdaptiveLocalRF = false;

    private String Filename;


    public class Model {
        @DefineFactor(onObservations = true)
        public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood =
                UnrootedTreeLikelihood.fromFastaFile(sequencesFile, selectedRateMtx).withExpFamMixture(ExpFamMixture.rateMtxModel(selectedRateMtx, isNormalized)).withTree(treeFile);

        @DefineFactor
        public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> priorOnParams =
                IIDRealVectorGenerativeFactor.iidNormalOn(likelihood.evolutionaryModel.rateMatrixMixture.parameters);
    }

    private final PrintWriter detailWriter = BriefIO.output(Results.getFileInResultFolder("experiment.details.txt"));

    public Model model;

    @Override
    public void run() {
        model = new Model();
        factory.addProcessor(this);
        factory.addProcessor(new LogDensityProcessor());

        long startTime = System.currentTimeMillis();

        if (useMH)
        {
            factory.excludeNodeMove(RealVectorOverRelaxedSlice.class);
            factory.excludeNodeMove(RealVectorAdaptiveMHProposal.class);

        }
        else factory.excludeNodeMove(RealVectorMHProposal.class);

        if (useGlobalRF) {
            factory.addNodeMove(ExpFamParameters.class, PhyloRFMove.class);
            factory.excludeNodeMove(RealVectorOverRelaxedSlice.class);
            factory.excludeNodeMove(RealVectorAdaptiveMHProposal.class);
            factory.excludeNodeMove(PhyloLocalRFMove.class);
            PhyloRFMove.usePegasusSolver = usePegasusSolver;
            PhyloRFMove.useSelfImplementedSolver = useSelfImplementedSolver;
            MultipleConvexCollisionSolver.useSelfImplementedSolver = useSelfImplementedSolver;

            if(usePegasusSolver && !useSelfImplementedSolver){
                PhyloRFMove.solverNames = SolverNames.Pegasus;
            }
            else if(!usePegasusSolver && useSelfImplementedSolver){
                PhyloRFMove.selfSolverNames = selfSolverNames;
            }
            else if(!useSelfImplementedSolver && !usePegasusSolver){
                PhyloRFMove.solverNames = SolverNames.Brent;
            }else{
                throw new RuntimeException("do not use Pegasus/Brent solvers and self implemented solvers at the same time");
            }

        }

        if (useLocalRF) {
            factory.addNodeMove(ExpFamParameters.class, PhyloLocalRFMove.class);

            if(useAdaptiveLocalRF){

                PhyloLocalRFMove.lowerbound = lowerbound;
                PhyloLocalRFMove.upperbound = upperbound;


            }else{

                if(maxTrajectoryLength == null){
                    throw new RuntimeException("The length of the trajectory should be provided if adaptation is not used");
                }
                PhyloLocalRFMove.maxTrajectoryLength = maxTrajectoryLength;
            }
            factory.excludeNodeMove(RealVectorOverRelaxedSlice.class);
            factory.excludeNodeMove(RealVectorAdaptiveMHProposal.class);
            factory.excludeNodeMove(PhyloRFMove.class);
        }

        else ;

        if (useHMC) ;
        else factory.excludeNodeMove(PhyloHMCMove.class);

        if(useAdaptiveHMC){
            factory.excludeNodeMove(RealVectorAdaptiveMHProposal.class);
            factory.excludeNodeMove(RealVectorOverRelaxedSlice.class);
            factory.excludeNodeMove(PhyloLocalRFMove.class);
            factory.excludeNodeMove(PhyloRFMove.class);
        }


        int nMovesRequested = (useGlobalRF ? 1 : 0) + (useHMC ? 1 : 0) + (useMH ? 1 : 0) + (useLocalRF ? 1 : 0) ;

        MCMCAlgorithm mcmc = factory.build(model, false);
        mcmc.options.nMCMCSweeps = nMCMCIterations;
        mcmc.options.burnIn=0;
        mcmc.options.thinningPeriod=1;

        if (mcmc.sampler.moves.size() != nMovesRequested)
            throw new RuntimeException("" + mcmc.sampler.moves.size() + "!=" + nMovesRequested);

        if (!useAdaptiveHMC && useHMC) {
            // find the hmc sampler
            PhyloHMCMove hmcMove = null;
            for (Move move : mcmc.sampler.moves)
                if (move instanceof PhyloHMCMove)
                    hmcMove = (PhyloHMCMove) move;
            // disable adaptivity by setting fixed L, epsilon values
            hmcMove.epsilon = this.epsilon;
            hmcMove.L = this.L;
        }

        if (useGlobalRF)
            for (Move move : mcmc.sampler.moves)
                if (move instanceof PhyloRFMove) {
                    ((PhyloRFMove) move).options = this.rfOptions;
                    ((PhyloRFMove) move).nItersPerPathAuxVar = this.nItersPerPathAuxVar;
                    ((PhyloRFMove) move).useSuperPosition = this.useSuperPosition;

                }

        if (useLocalRF)
            for (Move move : mcmc.sampler.moves)
                if (move instanceof PhyloRFMove) {
                    ((PhyloLocalRFMove) move).options = this.rfOptions;
                }

        if (useHMC)
            for (Move move : mcmc.sampler.moves)
                if (move instanceof PhyloHMCMove)
                    ((PhyloHMCMove) move).nItersPerPathAuxVar = this.nItersPerPathAuxVar;


        System.out.println(mcmc.model);
        System.out.println(mcmc.sampler);

        mcmc.run();

        logToFile("Total time in minutes: " + ((System.currentTimeMillis() - startTime)/60000.0));
        logToFile("Fixed Trajectory length:" + maxTrajectoryLength);
       // File newDirectory = new File(Results.getResultFolder().getParent() + "rep"+ rep+ "isExcludedHMCMove" + isExcluded + bandwidth+selectedRateMtx+"numSites"+numberOfSites+"Seed"+whichSeedUsed+ "epsilon"+PhyloHMCMove.epsilon+"L"+PhyloHMCMove.L);
        logToFile("Samplers:"+ mcmc.toString());

        if(useAdaptiveLocalRF){

            Filename = Results.getResultFolder().getParent() + "nIter"+nMCMCIterations+ "TrajLength"+ PhyloLocalRFMove.maxTrajectoryLength + "useLocal"+ useLocalRF +
                     "useAdaptiveLocal"+ useAdaptiveLocalRF + selectedRateMtx;
        }

        if(!useAdaptiveLocalRF&& useLocalRF){

            Filename = Results.getResultFolder().getParent() + "nIter"+nMCMCIterations+ "TrajLength"+ maxTrajectoryLength + "useLocal"+ useLocalRF +
                       selectedRateMtx;
        }
        if(useHMC||useAdaptiveHMC){

            Filename = Results.getResultFolder().getParent() + "nIter"+nMCMCIterations+ "epsilon"+ PhyloHMCMove.epsilon + "L"+ PhyloHMCMove.L+
                     selectedRateMtx;
        }

        if(useGlobalRF){
            Filename = Results.getResultFolder().getParent() + "nIter"+nMCMCIterations+ "useGlobal" + useGlobalRF + "usePegasusSolver" + usePegasusSolver + "useSelfImplementedSolver"+
                    useSelfImplementedSolver   +  "useSuperPosition"+ useSuperPosition + selectedRateMtx;

            if(usePegasusSolver && !useSelfImplementedSolver){
                Filename = Filename + "solverNames" + SolverNames.Pegasus;
            }
            else if(!usePegasusSolver && useSelfImplementedSolver){
                Filename = Filename + "solverNames" + selfSolverNames;
            }
            else if(!useSelfImplementedSolver && !usePegasusSolver){
                Filename = Filename + "solverNames" + SolverNames.Brent;
            }else{
                throw new RuntimeException("do not use Pegasus/Brent solvers and self implemented solvers at the same time");
            }

        }

        File newDirectory = new File(Filename);

        newDirectory.mkdir();
        try
        {
            FileUtils.copyDirectory(Results.getResultFolder(), newDirectory);
        } catch (IOException e)
        {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        Mains.instrumentedRun(args, new TestSimuDataLocalSampler());
    }

    Map<String, List<Double>> data = Maps.newLinkedHashMap();

    public void logToFile(String someline) {
        this.detailWriter.println(someline);
        this.detailWriter.flush();
    }


    @Override
    public void process(ProcessorContext context) {
    }
}
