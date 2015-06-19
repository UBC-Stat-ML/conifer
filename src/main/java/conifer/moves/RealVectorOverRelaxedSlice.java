package conifer.moves;

import bayonet.distributions.Normal;
import blang.factors.Factor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.NodeMove;
import blang.mcmc.SampledVariable;
import briefj.BriefIO;
import briefj.Indexer;
import briefj.opt.Option;
import briefj.run.Results;
import com.google.common.collect.Lists;
import conifer.ctmc.PathStatistics;
import conifer.ctmc.expfam.*;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import org.jblas.DoubleMatrix;
import org.apache.commons.math3.util.Pair;

import java.io.PrintWriter;
import java.util.List;
import java.util.Random;

/**
 * Created by crystal on 2015-06-16.
 * @author Tingting Zhao (zhaott0416@gmail.com)
 */
public class RealVectorOverRelaxedSlice extends NodeMove
{
    @SampledVariable
    ExpFamParameters parameters;

    @ConnectedFactor
    UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood;
    @ConnectedFactor
    IIDRealVectorGenerativeFactor<Normal.MeanVarianceParameterization> prior;

    private double [] savedValue = null;

    private final PrintWriter detailWriter = BriefIO.output(Results.getFileInResultFolder("SliceSampler.experiment.details.txt"));

    private final double SLICE_SIZE = 2;
    private final int MAX_SLICE_SIZE = 5; // max size of slice is MAX_SLICE_SIZE * SLICE_SIZE

    @Override
    public void execute(Random rand)
    {
        if (savedValue != null)
            throw new RuntimeException();
        if (prior.marginalDistributionParameters.mean.getValue() != 0.0)
            throw new RuntimeException();
        final double variance = prior.marginalDistributionParameters.variance.getValue();

        double [] initialPoint = parameters.getVector();

        savedValue = initialPoint.clone();

        for(int i=0; i<initialPoint.length; i++){

            double originalUnnormalizedLogLikelihood = likelihood.logDensity();
            double auxVariable = originalUnnormalizedLogLikelihood + Math.log(rand.nextDouble());

            double originalValue = savedValue[i];
            Pair<Double, Double> interval =steppingOutInterval(rand, auxVariable, originalValue, i);
            double newValue = shrinkingSampling(auxVariable, originalValue, interval, i);
            savedValue[i] = newValue;
            likelihood.evolutionaryModel.rateMatrixMixture.parameters.setVector(savedValue);
        }

        parameters.setVector(savedValue);
        savedValue = null;

    }



    private Pair<Double, Double> steppingOutInterval(Random rand, double auxVariable, double originalValue, int index)
    {
        // index is the location of originalValue in the original vector savedValue
        double L = originalValue - SLICE_SIZE * rand.nextDouble();
        double R = L + SLICE_SIZE;
        int J = (int) Math.floor(MAX_SLICE_SIZE * rand.nextDouble());
        int K = (MAX_SLICE_SIZE - 1) - J;

        double [] LVector= savedValue.clone();
        double [] RVector = savedValue.clone();
        LVector[index] = L;
        double leftLogUnnormalizedPotential = computeLogUnnormalizedPotentials(LVector);

        RVector[index] = R;
        double rightLogUnnormalizedPotential = computeLogUnnormalizedPotentials(RVector);

        while (J > 0 && (auxVariable < leftLogUnnormalizedPotential))
        {
            L -= SLICE_SIZE;
            J -= 1;
            leftLogUnnormalizedPotential = computeLogUnnormalizedPotentials(LVector);
        }

        while (K > 0 && (auxVariable < rightLogUnnormalizedPotential))
        {
            R += SLICE_SIZE;
            K -= 1;
            rightLogUnnormalizedPotential = computeLogUnnormalizedPotentials(RVector);
        }

        return new Pair<Double, Double>(L, R);
    }

    public double shrinkingSampling(double auxVariable, double originalValue, Pair<Double, Double> interval, int index)
    {
        double L = interval.getFirst();
        double R = interval.getSecond();
        double w = SLICE_SIZE;
        int a = MAX_SLICE_SIZE;


        if (R - L < 1.1 * w)
        {
            while(true)
            {
                double M = (L + R) / 2;
                double [] curVector = savedValue.clone();
                curVector[index] = M;

                if (a == 0 || auxVariable < computeLogUnnormalizedPotentials(curVector))
                {
                    break;
                }
                if (originalValue > M)
                {
                    L = M;
                }
                else
                {
                    R = M;
                }
                a -= 1;
                w *= 0.5;
            }
        }

        double Lhat = L;
        double Rhat = R;

        while(a > 0)
        {
            a -= 1;
            w *= 0.5;
            double [] LVector = savedValue.clone();
            double [] RVector = savedValue.clone();
            LVector[index] = Lhat+w;
            Lhat = (auxVariable >= computeLogUnnormalizedPotentials(LVector)) ? Lhat + w : Lhat;
            RVector[index] = Rhat-w;
            Rhat = (auxVariable >= computeLogUnnormalizedPotentials(RVector)) ? Rhat - w : Rhat;
        }

        double proposed = Lhat + Rhat - originalValue;


        double [] propVector = savedValue.clone();
        propVector[index] = proposed;
        if (proposed < L || proposed > R || auxVariable >= computeLogUnnormalizedPotentials(propVector))
        {
            proposed = originalValue;
        }

        return proposed;

    }


    public double computeLogUnnormalizedPotentials(double [] value)
    {
        likelihood.evolutionaryModel.rateMatrixMixture.parameters.setVector(value);
        return likelihood.logDensity();

    }



    public static List<CTMCState> states(int category, CTMCStateSpace space)
    {
        List<CTMCState> result = Lists.newArrayList();

        Object partitionId = space.currentPartition;
        Indexer<?> latentIndexer = space.latentIndexer;

        for (int i = 0; i < latentIndexer.size(); i++)
            result.add(new CTMCState(category, latentIndexer.i2o(i), partitionId));

        return result;
    }

    public void logToFile(String someline) {
        this.detailWriter.println(someline);
        this.detailWriter.flush();
    }


    @Override
    public String toString()
    {
        return "RealVectorOverRelaxedSlice [parameters=" + parameters  + "]";
    }





}
