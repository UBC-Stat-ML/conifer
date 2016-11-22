package conifer.rejfreemodels.phylo;


import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import blang.ProbabilityModel;
import briefj.OutputManager;
import briefj.run.Results;
import com.google.common.base.Stopwatch;
import conifer.global.SelfImplementedSolver;
import conifer.global.SelfImplementedSolverNames;
import conifer.global.SolverNames;
import org.jblas.DoubleMatrix;

import conifer.rejfreeutil.RFSamplerOptions;
import conifer.global.GlobalRFSampler;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.NodeMove;
import blang.mcmc.SampledVariable;
import conifer.ctmc.PathStatistics;
import conifer.ctmc.expfam.CTMCExpFam;
import conifer.ctmc.expfam.CTMCState;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.ExpFamParameters;
import conifer.ctmc.expfam.ExpectedStatistics;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.moves.PhyloHMCMove;


public class PhyloRFMove extends NodeMove {
    @SampledVariable
    ExpFamParameters parameters;

    @ConnectedFactor
    UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood;
    @ConnectedFactor
    IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior;

    public RFSamplerOptions options = new RFSamplerOptions();

    public int nItersPerPathAuxVar = 1000;

    public ExpectedCompleteReversibleModel modelSpec;

    public static boolean useSuperPosition = true;

    public static boolean usePegasusSolver = true;

    public static boolean useSelfImplementedSolver = false;

    public static SolverNames solverNames = SolverNames.Pegasus;

    public static SelfImplementedSolverNames selfSolverNames = SelfImplementedSolverNames.ImprovedPegasus;


    public Stopwatch watch = null;

    public OutputManager output = Results.getGlobalOutputManager();

    @Override
    public void execute(Random rand) {
        if (prior.marginalDistributionParameters.mean.getValue() != 0.0)
            throw new RuntimeException();
        final double variance = prior.marginalDistributionParameters.variance.getValue();

        List<PathStatistics> pathStatistics = likelihood.evolutionaryModel.samplePosteriorPaths(rand, likelihood.observations, likelihood.tree);

        ExpectedStatistics<CTMCState> convertedStat = PhyloHMCMove.convert(pathStatistics, parameters, likelihood);
        CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective objective = parameters.globalExponentialFamily.getExpectedCompleteReversibleObjective(1.0 / variance, convertedStat);
        modelSpec = new ExpectedCompleteReversibleModel(parameters, objective, parameters.globalExponentialFamily);


        double[] initialPoint = parameters.getVector();


        GlobalRFSampler sampler;

        if (initialized) {
            sampler = new GlobalRFSampler(objective, new DoubleMatrix(initialPoint), options, new ProbabilityModel(modelSpec));
        } else {
            System.out.println("Initializing RF sampler");

            if(!useSelfImplementedSolver){
                sampler = GlobalRFSampler.initializeRFWithLBFGS(objective, options, new ProbabilityModel(modelSpec), solverNames);
            }else{

                sampler = GlobalRFSampler.initializeRFWithLBFGS(objective, options, new ProbabilityModel(modelSpec), selfSolverNames);
            }

//            if(usePegasusSolver && !useSelfImplementedSolver){
//                sampler = GlobalRFSampler.initializeRFWithLBFGS(objective, options, new ProbabilityModel(modelSpec));
//            }
//
//            else if( useSelfImplementedSolver && !usePegasusSolver){
//
//                sampler = GlobalRFSampler.initializeRFWithLBFGS(objective, options, new ProbabilityModel(modelSpec), selfSolverNames);
//            }
//
//            else if(!usePegasusSolver && !useSelfImplementedSolver){
//                sampler = GlobalRFSampler.initializeRFWithLBFGS(objective, options, new ProbabilityModel(modelSpec), SolverNames.Brent);
//            }
//
//            else{
//                throw new RuntimeException("We do not use Pegasus/Brent together with self implemented solvers simutaneously");
//
//            }
            initialized = true;
        }

        watch = Stopwatch.createStarted();

        if(!useSuperPosition){
            sampler.iterate(rand, nItersPerPathAuxVar);
        }else{
            sampler.iterateWithSuperPosition(rand, nItersPerPathAuxVar, true);
        }

        watch.stop();
        long elapsed = watch.elapsed(TimeUnit.MILLISECONDS);


        double[] newPoint = sampler.getCurrentPosition().data;

        parameters.setVector(newPoint);

        if(useSuperPosition){
            output.printWrite("general-sampler-diagnostic",
                    "wallClockTimeMilliSeconds", elapsed,
                    "totalTime", sampler.getTotalTimeWithSuperposition(),
                    "nBounces", sampler.getNBounces(),
                    "nVirtualCollisions", sampler.getNVirtualCollisions(),
                    "collisionCalculationTime", sampler.getCollisionCalculationTime(),
                    "candidateCollisionCalculationTime", sampler.getCandidateCollisionCalculationTime(),
                    "gradientCalculationTime", sampler.getGradientCalculationTime(),
                    "nRefreshment", sampler.getNRefreshments());
        }else{
            output.printWrite("general-sampler-diagnostic",
                    "wallClockTimeMilliSeconds", elapsed,
                    "averageTime", sampler.getAverageTimeWithoutSuperposition(),
                    "nBounces", sampler.getNBounces(),
                    "collisionCalculationTime", sampler.getCollisionCalculationTime(),
                    "gradientCalculationTime", sampler.getGradientCalculationTime(),
                    "refreshmentTime", sampler.getRefreshmentTime(),
                    "nRefreshment", sampler.getNRefreshments());
        }

        output.flush();
    }

    private boolean initialized = false;
}
