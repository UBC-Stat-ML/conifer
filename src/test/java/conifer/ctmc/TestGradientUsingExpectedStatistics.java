package conifer.ctmc;

import bayonet.distributions.Normal;
import bayonet.math.EJMLUtils;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import briefj.run.Results;
import com.beust.jcommander.internal.Lists;
import com.google.gson.Gson;
import conifer.SimplePhyloSimulatorFixedTreeRateMtx;
import conifer.ctmc.expfam.*;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.io.PhylogeneticObservationFactory;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.moves.PhyloHMCMove;
import org.ejml.simple.SimpleMatrix;
import org.jblas.DoubleMatrix;
import org.junit.Test;


import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import briefj.BriefCollections;
import briefj.BriefIO;
import briefj.BriefLog;

/**
 * Test the agreement between pure forward sampler vs.
 * forward+end point sampler.
 *
 * @author Tingting Zhao (zhaott0416@gmail.com)
 *
 */
public class TestGradientUsingExpectedStatistics
{


    private final PrintWriter detailWriter = BriefIO.output(Results.getFileInResultFolder("experiment.details.txt"));

    public static void main(String [] args)
    {
        new TestGradientUsingExpectedStatistics().test();
    }

    private final double Threshold = 0.5;

    @Test
    public void test()
    {
//        File inputFile = new File("src/main/resources/conifer/sampleInput/FES_4.fasta");
//        File treeFile = new File("src/main/resources/conifer/sampleInput/FES.ape.4.nwk");
//        RateMtxNames selectedRateMtx = RateMtxNames.DNAGTR;
//        CTMCExpFam<CTMCState>.ExpectedReversibleObjectiveUpdateExpectedStat noAuxObjective =null;
//        class Model {
//
//            @DefineFactor(onObservations = true)
//            public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood1 =
//                    UnrootedTreeLikelihood
//                            .fromFastaFile(inputFile, selectedRateMtx)
//                            .withExpFamMixture(ExpFamMixture.rateMtxModel(selectedRateMtx, false))
//                            .withTree(treeFile);
//
//            @DefineFactor
//            public final IIDRealVectorGenerativeFactor<Normal.MeanVarianceParameterization> prior =
//                    IIDRealVectorGenerativeFactor
//                            .iidNormalOn(likelihood1.evolutionaryModel.rateMatrixMixture.parameters);
//
//        }
//
//        Model model = new Model();
//
//        UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = model.likelihood1;
//        ExpFamParameters parameters = likelihood.evolutionaryModel.rateMatrixMixture.parameters;
//
//
//        int numOfWeights = 5;
//        for(int j=0; j< numOfWeights; j++){
//            double [] newWeight = new double [parameters.getVector().length];
//            for(int i=0; i< parameters.getVector().length;i++){
//                newWeight[i] = new Random(i+j).nextGaussian();
//            }
//
//            parameters.setVector(newWeight);
//
//            System.out.println(Arrays.toString(likelihood.evolutionaryModel.rateMatrixMixture.parameters.getVector()));
//
//            ExpectedStatistics<CTMCState> expectedStatistics = likelihood.evolutionaryModel.getTotalExpectedStatistics(likelihood.observations, likelihood.tree, parameters.globalExponentialFamily);
//            noAuxObjective = parameters.globalExponentialFamily.getExpectedReversibleObjectiveUpdateExpectedStat(1.0, likelihood.observations,
//                    likelihood.tree, likelihood.evolutionaryModel, expectedStatistics, parameters);
//            double [] gradientFromExpectedStat = noAuxObjective.derivativeAt(parameters.getVector());
//
//            // next I get the averaged gradient using auxiliary variables and average the gradients obtained from sampled values of the sufficient statistics
//            List<List<PathStatistics>> allPathStatistics = Lists.newArrayList();
//            List<ExpectedStatistics<CTMCState>> convertedStatLists = Lists.newArrayList();
//            List<CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective> auxObjectiveLists = Lists.newArrayList();
//
//            int nRep =5000;
//            DoubleMatrix allInit = new DoubleMatrix(1, likelihood.evolutionaryModel.rateMatrixMixture.getRateMatrix(0).getRateMatrix().length);
//            DoubleMatrix allHoldingtime = new DoubleMatrix(1, likelihood.evolutionaryModel.rateMatrixMixture.getRateMatrix(0).getRateMatrix().length);
//            DoubleMatrix allNTrans = new DoubleMatrix(likelihood.evolutionaryModel.rateMatrixMixture.getRateMatrix(0).getRateMatrix().length,
//                    (likelihood.evolutionaryModel.rateMatrixMixture.getRateMatrix(0).getRateMatrix().length-1));
//
//
//            for(int i=0; i< nRep; i++){
//                allPathStatistics.add(i, likelihood.evolutionaryModel.samplePosteriorPaths(new Random(i), likelihood.observations, likelihood.tree));
//                convertedStatLists.add(i, PhyloHMCMove.convert(allPathStatistics.get(i), parameters, likelihood));
//                // sum over all the sufficient statistics across all iterations and outside this for loop, we will get an average value
//                allInit.addi(new DoubleMatrix(convertedStatLists.get(i).nInit).transpose());
//                allHoldingtime.addi(new DoubleMatrix((convertedStatLists.get(i).holdTimes)).transpose());
//                allNTrans.addi(new DoubleMatrix(convertedStatLists.get(i).nTrans));
//                //auxObjectiveLists.add(i, parameters.globalExponentialFamily.getExpectedCompleteReversibleObjective(1.0, convertedStatLists.get(i)));
//            }
//
//            // get averaged sufficient statistics of NInit, holdingTime and nTrans across all nRep iterations
//            double [] avgInit = allInit.div(nRep).toArray();
//            double [] avgHoldingTime = allHoldingtime.div(nRep).toArray();
//            double [][] avgNTrans = allNTrans.div(nRep).toArray2();
//
//            ExpectedStatistics<CTMCState> avgExpectedStat = convertedStatLists.get(0);
//            avgExpectedStat.nInit = avgInit;
//            avgExpectedStat.holdTimes=avgHoldingTime;
//            avgExpectedStat.nTrans=avgNTrans;
//
//            CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective auxObjective =  parameters.globalExponentialFamily.getExpectedCompleteReversibleObjective(1.0, avgExpectedStat);
//            double [] avgGradient = auxObjective.derivativeAt(parameters.getVector());
//            System.out.println(Arrays.toString(avgGradient));
//
////            double [][] allGradient = new double [nRep][parameters.getVector().length];
////            for(int i=0; i< nRep; i++){
////                allGradient[i] = auxObjectiveLists.get(i).derivativeAt(parameters.getVector());
////            }
////            DoubleMatrix allGradientMtx = new DoubleMatrix(allGradient);
////            double []  averagedGradient = allGradientMtx.columnMeans().toArray();
////            System.out.println(Arrays.toString(averagedGradient));
//            System.out.println(Arrays.toString(gradientFromExpectedStat));
////            EJMLUtils.checkIsClose(new SimpleMatrix(allGradientMtx.columnMeans().transpose().toArray2()), new SimpleMatrix(new DoubleMatrix(gradientFromExpectedStat).toArray2()), 0.5);
//            EJMLUtils.checkIsClose(new SimpleMatrix(new DoubleMatrix(avgGradient).toArray2()), new SimpleMatrix(new DoubleMatrix(gradientFromExpectedStat).toArray2()), Threshold);
//            System.out.println(j);
//        }



    }


}
