package conifer.rejfreemodels.phylo;


import java.util.List;
import java.util.Random;

import blang.ProbabilityModel;
import briefj.opt.OptionSet;
import conifer.local.LocalRFRunner;
import conifer.local.LocalRFRunnerOptions;
import conifer.local.LocalRFSampler;
import conifer.rejfreemodels.normal.ComparisonUtils;
import conifer.rejfreeprocessors.MomentRayProcessor;
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


    @Override
    public void execute(Random rand) {
        if (prior.marginalDistributionParameters.mean.getValue() != 0.0)
            throw new RuntimeException();
        final double variance = prior.marginalDistributionParameters.variance.getValue();

        List<PathStatistics> pathStatistics = likelihood.evolutionaryModel.samplePosteriorPaths(rand, likelihood.observations, likelihood.tree);

        ExpectedStatistics<CTMCState> convertedStat = PhyloHMCMove.convert(pathStatistics, parameters, likelihood);
        CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective objective = parameters.globalExponentialFamily.getExpectedCompleteReversibleObjective(1.0 / variance, convertedStat);
        modelSpec = new ExpectedCompleteReversibleModel(parameters, objective, parameters.globalExponentialFamily);

        localRFRunnerOption.maxSteps = Integer.MAX_VALUE;
        localRFRunnerOption.maxTrajectoryLength = 0.01;
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
}
