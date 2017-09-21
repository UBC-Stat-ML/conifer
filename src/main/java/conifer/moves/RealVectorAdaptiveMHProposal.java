package conifer.moves;

import blang.MCMCAlgorithm;
import blang.factors.Factor;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.MHProposalDistribution;
import blang.mcmc.SampledVariable;
import blang.processing.ProcessorContext;
import blang.variables.RealVectorInterface;
import org.jblas.Decompose;
import org.jblas.DoubleMatrix;

import java.util.List;
import java.util.Random;

// TODO: this should go in blang!

/**
 * @author Tingting Zhao (zhaott0416@gmail.com)
 * This class implements section 2 "Adaptive Metropolis" of paper "Examples of Adaptive MCMC" by Gareth O Roberts and Jeffrey S. Rosenthal
 * This is an adaptive MCMC algorithm since after the MCMC iteration is bigger than 2*dim(dimension of the parameters), we propose parameters from
 * a mixture of Normal kernels with weight beta suggested as 0.05 in the paper, the first kernel is N(x, 2.38^2/dim * sigma_n), where sigma_n is
 * the empirical sample covariance matrix for all previous MCMC iterations, x is current value of the parameter. The second kernel is N(x, (0.1)^2I_d/dim)
 * I_d is the identity matrix of dimension m. In the implementation below, we use the recursion formula to calculate the mean of all previous samples and
 * the covariance matrix from iteration i to iteration (i+1) since we only do not keep all previous values of the parameters in our implementation.
 */

public class RealVectorAdaptiveMHProposal implements MHProposalDistribution
{
    public static double bandWidth1 = 0.01;
    public static double bandWidth2 = 2.38*2.38;

    // this counter records the number of times that the propose method is called, it reflects which
    // MCMC iteration the sampler is at.
    public static int counter =0;

    @SampledVariable RealVectorInterface variable;

    @ConnectedFactor List<Factor> connectedFactors;

    private double [] savedValue = null;


    // meanVector is the mean of the variable across all previous iterations
    public static double [] meanVector = null;

    // The proposal is Qn, which is a mixture of two kernels. Qn=(1-beta)* N(x, (2.38)^2* sampleCov/ d) + beta* N(x, (0.1)^2* I_d/d), d is the dimension and sampleCov is the
    // empirical covariance matrix for all previous iterations. This proposal is a mixture of two kernels. If n<=2d, use the second
    // kernel, otherwise use this mixture of kernels.
    public static double [][] sampleCov = null;
    public static double beta=0.05;

    @Override
    public Proposal propose(Random rand)
    {
        //System.out.println("Computing RealVectorHMProposal");
        if (savedValue != null)
            throw new RuntimeException();
        double [] variableArray = variable.getVector();
        savedValue = variableArray.clone();

        //propose new value based on a multivariate Normal distribution with mean and covariance
        int dim = variable.getDim();

        DoubleMatrix variableArrayMtx = new DoubleMatrix(variableArray);

        counter++;

        if(counter<=2*dim){
            variableArray =  generateMultiNormal(rand, variableArray, DoubleMatrix.eye(dim).toArray2(), Math.sqrt(bandWidth1/dim));

        }else{

            DoubleMatrix covariance = new DoubleMatrix(sampleCov);
            double tmp = rand.nextDouble();
            if(tmp<0.05){
                variableArray =  generateMultiNormal(rand, variableArray, DoubleMatrix.eye(dim).toArray2(), Math.sqrt(bandWidth1/dim));

            }else{
                variableArray = generateMultiNormal(rand, variableArray, sampleCov, Math.sqrt(bandWidth2 / dim));

            }
        }
        variable.setVector(variableArray);

        return new ProposalRealization();
    }



    public static double [] generateMultiNormal(Random rand, double [] mean,  double [][]cov, double bandWidth){

        // bandwidth is used when we generate the next state from normal x' = x+ bandwidth * U* standardNormal, U is the cholesky
        // decomposition of the covaraince matrix
        int dim = mean.length;
        // generate multivariate standard normal first
        double [] standardMultiNormal = new double [dim];
        for(int i=0; i< dim; i++){
            standardMultiNormal[i] = rand.nextGaussian();
        }

        DoubleMatrix standardMultiNormalMtx = new DoubleMatrix(standardMultiNormal);
        DoubleMatrix meanMtx = new DoubleMatrix(mean);
        DoubleMatrix covMtx = new DoubleMatrix(cov);

        // get cholesky decomposition of the covariance matrix
        DoubleMatrix scaleMtx = Decompose.cholesky(covMtx).transpose(); // scaleMtx is a lower triangle matrix
        double [] result;
        result = meanMtx.add(scaleMtx.mmul(standardMultiNormalMtx).mul(bandWidth)).toArray();
        return result;
    }


    public static void updateMeanAndCov(double curIter, double [] curVector){

        // Need to update the meanVector and sampleCov based on the recursion
        // I need to pass the iteration number of MCMC to the proposal
        DoubleMatrix curVectorMtx = new DoubleMatrix(curVector);
        double [] oldMeanVector = meanVector.clone();
        DoubleMatrix oldMeanVectorMtx = new DoubleMatrix(oldMeanVector);
        meanVector = oldMeanVectorMtx.mul(curIter-1).add(curVectorMtx).div(curIter).toArray();
        DoubleMatrix meanVectorMtx = new DoubleMatrix(meanVector);

        if(curIter==1){
            sampleCov = DoubleMatrix.eye(curVector.length).toArray2();
        }else{

            //To update the sampleCov, there are mainly four components,
            DoubleMatrix tmp1, tmp2, tmp3, tmp4;
            tmp1 = new DoubleMatrix(sampleCov).mul((curIter - 1) / curIter);
            tmp2 = oldMeanVectorMtx.mmul(oldMeanVectorMtx.transpose());
            tmp3 = meanVectorMtx.mmul(meanVectorMtx.transpose()).mul((-1)*(curIter)/(curIter-1));
            tmp4 = curVectorMtx.mmul(curVectorMtx.transpose()).div(curIter-1);
            sampleCov = tmp1.add(tmp2).add(tmp3).add(tmp4).toArray2();

        }
    }



    private class ProposalRealization implements Proposal
    {
        private double curIter;

        @Override
        public double logProposalRatio()
        {
            return 0;
        }

        @Override
        public void acceptReject(boolean accept) {
            curIter = (double) counter;

            if (!accept) {
                double[] variableArray = variable.getVector();
                for (int i = 0; i < variableArray.length; i++)
                    variableArray[i] = savedValue[i];
                variable.setVector(variableArray);

            }

            updateMeanAndCov(curIter, variable.getVector());

            savedValue = null;
        }

    }




}
