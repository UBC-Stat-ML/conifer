package conifer.rejfreemodels.phylo;


import java.io.PrintWriter;
import java.util.List;
import java.util.Random;

import blang.ProbabilityModel;
import briefj.BriefIO;
import briefj.opt.OptionSet;
import briefj.run.Results;
import conifer.local.ARFS;
import conifer.local.LocalRFRunner;
import conifer.local.LocalRFRunnerOptions;
import conifer.local.LocalRFSampler;
import conifer.rejfreemodels.normal.ComparisonUtils;
import conifer.rejfreeprocessors.MomentRayProcessor;
import hmc.AHMC;
import hmc.DataStruct;
import hmc.HMC;
import org.apache.commons.lang3.time.StopWatch;
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


public class PhyloLocalRFMove extends NodeMove {
    @SampledVariable
    ExpFamParameters parameters;

    @ConnectedFactor
    UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood;
    @ConnectedFactor
    IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior;

    public RFSamplerOptions options = new RFSamplerOptions();

    public LocalRFRunnerOptions localRFRunnerOption = new LocalRFRunnerOptions();

    public LocalRFSampler sampler;

    public ExpectedCompleteReversibleModel modelSpec;

    public static Double maxTrajectoryLength = null;

    private final PrintWriter detailWriter = BriefIO.output(Results.getFileInResultFolder("RejectionFree.experiment.details.txt"));

    public static boolean useAdaptiveTrajectoryLength = true;

    public static Double lowerbound = null;
    public static Double upperbound = null;


    @Override
    public void execute(Random rand) {
        if (prior.marginalDistributionParameters.mean.getValue() != 0.0)
            throw new RuntimeException();
        final double variance = prior.marginalDistributionParameters.variance.getValue();
        double [] initialPoint = parameters.getVector();

        List<PathStatistics> pathStatistics = likelihood.evolutionaryModel.samplePosteriorPaths(rand, likelihood.observations, likelihood.tree);

        ExpectedStatistics<CTMCState> convertedStat = PhyloHMCMove.convert(pathStatistics, parameters, likelihood);
        CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective objective = parameters.globalExponentialFamily.getExpectedCompleteReversibleObjective(1.0 / variance, convertedStat);
        modelSpec = new ExpectedCompleteReversibleModel(parameters, objective, parameters.globalExponentialFamily);

        localRFRunnerOption.maxSteps = Integer.MAX_VALUE;
        if(!hyperParametersInitialized())
        {
            StopWatch watch = new StopWatch();
            watch.start();
            ARFS arfs  = null;

            if(useAdaptiveTrajectoryLength){
                arfs = ARFS.initializeARFSWithLBFGS(10000, 1000, objective, objective, initialPoint.length,modelSpec, lowerbound, upperbound);
            }

            double [] newPoint = new double [parameters.getVector().length];
            newPoint = arfs.sample(rand).data;
            maxTrajectoryLength = arfs.getFixedTrajectoryLength();
            logToFile("time to get adaptive trajectory length of local rejection free sampler:" + watch.getTime()/1000);
            logToFile("optimal trajectory length" +""+ maxTrajectoryLength);

        }

        localRFRunnerOption.maxTrajectoryLength = maxTrajectoryLength;
        LocalRFRunner rfRunner = new LocalRFRunner(localRFRunnerOption);

        rfRunner.init(modelSpec);
        rfRunner.addMomentRayProcessor();
        rfRunner.run();

        int nVariables = rfRunner.model.getLatentVariables().size();
        double [] newPoints = new double[nVariables];
        for(int i=0; i < nVariables; i++){
            newPoints[i] = modelSpec.variables.list.get(i).getValue();
        }

        parameters.setVector(newPoints);

    }

    private boolean hyperParametersInitialized()
    {
        return (maxTrajectoryLength != null);
    }

    public void logToFile(String someline) {
        this.detailWriter.println(someline);
        this.detailWriter.flush();
    }
}
