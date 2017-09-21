package conifer.moves;

import bayonet.distributions.Normal;
import bayonet.math.EJMLUtils;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import briefj.BriefIO;
import briefj.run.Results;
import com.beust.jcommander.internal.Lists;
import conifer.ctmc.PathStatistics;
import conifer.ctmc.expfam.*;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.RandomMatrices;
import org.ejml.simple.SimpleMatrix;
import org.jblas.DoubleMatrix;
import org.junit.Test;

import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Test the agreement between pure forward sampler vs.
 * forward+end point sampler.
 *
 * @author Tingting Zhao (zhaott0416@gmail.com)
 *
 */
public class TestRealVectorAdaptiveMHProposal
{
    @Test
    public void test()
    {
        int dim =8;
        // generate a matrix of 100* 8, 100 observations of dimension 8;
        // Test generateMultiNormal and updateMeanAndCov method of this class

        //Test generateMultiNormal method, generate data first under an arbitrary covariance matrix
        Random rand = new Random(1);
        double [][] cov = new double[dim][dim];
        DenseMatrix64F covariance = RandomMatrices.createSymmPosDef(dim, rand);
        for(int i=0; i<dim; i++){
            for(int j=0; j< dim; j++){
                cov[i][j] = covariance.get(i, j);

            }
        }

        double [] mean = new double[dim];
        for(int i=0; i< dim; i++){
            mean[i] = rand.nextGaussian();
        }
        System.out.println(Arrays.toString(mean));

        int nSamples = 50000;
        double [][] samples = new double [nSamples][dim];

        for(int j=0; j< nSamples; j++){

            samples[j] = RealVectorAdaptiveMHProposal.generateMultiNormal(rand, mean, cov, 1);

        }

        // obtain sample mean estimate
        double [] meanEst = new double [dim];
        meanEst = new DoubleMatrix(samples).columnMeans().toArray();
        System.out.println(Arrays.toString(meanEst));

        // test if the true mean and sample mean are close enough
        EJMLUtils.checkIsClose(new SimpleMatrix(new DoubleMatrix(mean).toArray2()), new SimpleMatrix(new DoubleMatrix(meanEst).toArray2()), 0.05);

        // calculate the sample covariance matrix to see if it is close to the true covariance matrix cov
        double [][] sampleCov = new double [dim][dim];
        DoubleMatrix sampleCovMtx = new DoubleMatrix(sampleCov);
        DoubleMatrix sampleMtx = new DoubleMatrix(samples);
        DoubleMatrix meanEstMtx = new DoubleMatrix(meanEst);
        DoubleMatrix tmp;
        DoubleMatrix tmpCovEachObs;
        for(int i =0; i< nSamples; i++){
            tmp = sampleMtx.getRow(i).sub(meanEstMtx);
            tmpCovEachObs = tmp.transpose().mmul(tmp).div(nSamples);
            sampleCovMtx = sampleCovMtx.add(tmpCovEachObs);
        }


        sampleCov = sampleCovMtx.toArray2();
        EJMLUtils.checkIsClose(new SimpleMatrix(cov), new SimpleMatrix(sampleCov), 0.02);

        System.out.println(Arrays.deepToString(cov));
        System.out.println(Arrays.deepToString(sampleCov));


        // check the updateMeanAndCov() method, this is to check whether the sample mean and sample covariance matrix are the same
        // either we calculate them directly or with the recursion formula used in method RealVectorAdpativeMHProposal.updateMeanAndCov()
        RealVectorAdaptiveMHProposal.meanVector = samples[0];
        RealVectorAdaptiveMHProposal.sampleCov = DoubleMatrix.eye(dim).toArray2();

        for(int j= 1; j< nSamples; j++){

            RealVectorAdaptiveMHProposal.updateMeanAndCov(j, samples[j]);
        }

        System.out.println(Arrays.deepToString(RealVectorAdaptiveMHProposal.sampleCov));

        EJMLUtils.checkIsClose(new SimpleMatrix(sampleCov), new SimpleMatrix(RealVectorAdaptiveMHProposal.sampleCov), 0.001);

        double [][] meanEst2d = new double [1][meanEst.length];
        meanEst2d[0] = meanEst;
        double [][] meanVector = new double [1][RealVectorAdaptiveMHProposal.meanVector.length];
        meanVector[0] = RealVectorAdaptiveMHProposal.meanVector;

        System.out.println(Arrays.toString(meanEst));
        System.out.println(Arrays.toString(RealVectorAdaptiveMHProposal.meanVector));

        EJMLUtils.checkIsClose(new SimpleMatrix(meanEst2d), new SimpleMatrix(meanVector), 0.001);



    }


}
