package conifer.local;

/**
 * Created by crystal on 2016-10-18.
 */
import java.util.Random;

import conifer.rejfreemodels.phylo.ExpectedCompleteReversibleModel;
import conifer.rejfreeutil.RFSamplerOptions;
import hmc.AHMC;
import hmc.DataStruct;
import hmc.HMC;
import hmc.HMC;
import org.jblas.DoubleMatrix;
import org.jblas.Solve;
import utils.MultiVariateObj;
import utils.Objective;

import bayonet.opt.LBFGSMinimizer;
import bo.BayesOpt;
import bo.kernel.CovModel;
import bo.kernel.CovSEARD;
/**
 * Created by crystal on 2016-10-18.
 */
public class ARFS {
    private int burnIn = 0;
    private int numIterations = 0;
    private MultiVariateObj gradient = null;
    private Objective fun = null;
    private int D;
    private int sizeAdapt = 10;
    private DoubleMatrix bound;
    private double noise = 0.1;
    private DoubleMatrix hyper =
            new DoubleMatrix(new double[] {0.2, 5});
    private DoubleMatrix initPt;
    private BayesOpt bo;

    public DoubleMatrix samples = null;
    public boolean adjustReward = true;
    private double fixedTrajectoryLength;

    private RFSamplerOptions options = new RFSamplerOptions();

    private LocalRFRunnerOptions localRFRunnerOption = new LocalRFRunnerOptions();

    private LocalRFSampler sampler;

    private ExpectedCompleteReversibleModel modelSpec;

    public static Double lowerbound = null;
    public static Double upperbound = null;

    public static ARFS initializeARFSWithLBFGS(int numIterations, int burnIn, MultiVariateObj gradient,
                                               Objective func, int D, ExpectedCompleteReversibleModel modelSpec, double lowerbound, double upperbound) {
        return new ARFS(numIterations, burnIn, defaultTrajectoryLengthBounds(lowerbound, upperbound), gradient, func, null, D, modelSpec);
    }

    public ARFS(int numIterations, int burnIn, MultiVariateObj gradient,
                Objective func, double [] initialPoint, ExpectedCompleteReversibleModel modelSpec, double lowerbound, double upperbound) {
        this(numIterations, burnIn, defaultTrajectoryLengthBounds(lowerbound, upperbound), gradient, func, initialPoint, initialPoint.length, modelSpec);
    }

    public ARFS(int numIterations, int burnIn, DoubleMatrix bound,
                MultiVariateObj gradient, Objective func, double [] initialPoint, int D, ExpectedCompleteReversibleModel modelSpec) {
        this.burnIn = burnIn;
        this.numIterations = numIterations;
        this.sizeAdapt = (int)Math.floor(this.burnIn/100.0);
        this.bound = bound;
        this.D = D;
        this.gradient = gradient;
        this.fun = func;
        if (initialPoint == null)
            initialPoint = initPoint_pureJava();
        this.initPt = new DoubleMatrix(initialPoint).transpose();
        CovModel cov = new CovSEARD(hyper);

        this.modelSpec = modelSpec;
        this.bo = new BayesOpt(0.0, convert(initFixedTrajectoryLength()), cov, bound, noise);
        this.bo.setUseScale(true);
        this.bo.setSooIter(200);
    }



    /**
     * Call this before sample() to request all samples
     * be kept in a matrix.
     */
    public void keepSamples()
    {
        this.samples = DoubleMatrix.zeros(this.numIterations-this.burnIn, this.D);
    }

    private double [] initPoint_pureJava() {

        LBFGSMinimizer minimizer = new LBFGSMinimizer();

        bayonet.opt.DifferentiableFunction f = new bayonet.opt.DifferentiableFunction() {

            @Override
            public int dimension()
            {
                return ARFS.this.D;
            }

            @Override
            public double valueAt(double[] x)
            {
                return ARFS.this.fun.functionValue(new DoubleMatrix(x));
            }

            @Override
            public double[] derivativeAt(double[] x)
            {
                return ARFS.this.gradient.mFunctionValue(new DoubleMatrix(x)).data;
            }

        };

        return minimizer.minimize(f, new double[this.D], 1e-5);
    }

    private DoubleMatrix convert(double fixedTrajectoryLength) {
        DoubleMatrix param = new DoubleMatrix(new double[]{Math.log(fixedTrajectoryLength)});
        DoubleMatrix ptEval = param.sub(this.bound.getColumn(0))
                .div(this.bound.getColumn(1).sub(
                        this.bound.getColumn(0))).transpose();
        return ptEval;
    }

    private double initFixedTrajectoryLength() {
        return Math.exp(this.bound.rowMeans().toArray()[0]);
    }


    public DoubleMatrix sample(Random rand) {
        fixedTrajectoryLength = initFixedTrajectoryLength();
        int numAdapt = 0; double reward = 0;
        DoubleMatrix ptEval = convert(fixedTrajectoryLength);
        DoubleMatrix sample = this.initPt.transpose();

        for (int ii = 0; ii < this.numIterations; ii++) {
            if (ii % this.sizeAdapt == 0) {
                if (this.adjustReward) {reward = reward / fixedTrajectoryLength;}
                numAdapt = numAdapt + 1;
                if (ii >= this.sizeAdapt) {this.bo.updateModel(ptEval, reward);}
                double rate = anneal(ii);

                if (rand.nextDouble() < rate) {
                    DoubleMatrix nextPt = this.bo.maximizeAcq(rate).
                            mul(this.bound.getColumn(1).
                                    sub(this.bound.getColumn(0))).
                            add(this.bound.getColumn(0));
                    fixedTrajectoryLength = Math.exp(nextPt.toArray()[0]);
                    ptEval = convert(fixedTrajectoryLength);
                }
                System.out.format("Iter  %3d fixedtrajectoryLength: %f " +
                                "reward: %f prob: %f\n",
                        numAdapt,fixedTrajectoryLength, reward, rate);
                reward = 0;
            }

            //change this part to local rejection free sampler
            //remove result.mr
            localRFRunnerOption.maxSteps = Integer.MAX_VALUE;
            localRFRunnerOption.maxTrajectoryLength = fixedTrajectoryLength;
            LocalRFRunner rfRunner = new LocalRFRunner(localRFRunnerOption);

            rfRunner.init(modelSpec);
            rfRunner.addMomentRayProcessor();
            rfRunner.run();

            int nVariables = rfRunner.model.getLatentVariables().size();
            double [] newPoints = new double[nVariables];
            for(int i=0; i < nVariables; i++){
                newPoints[i] = modelSpec.variables.list.get(i).getValue();
            }

            ARFSDataStruct result = new ARFSDataStruct(new DoubleMatrix(newPoints),
                    sample, fixedTrajectoryLength);

            sample = result.next_q;
            reward = reward +  Math.pow(result.next_q.sub(result.q).norm2(), 2);

            if (ii >= this.burnIn) {
                if (this.samples != null)
                    this.samples.putRow(ii-this.burnIn, sample.transpose());
            }
        }
        return sample;
    }

    public double anneal(int ii) {
        double anneal = Math.pow(Math.max( (ii - this.burnIn) /
                ((double) this.sizeAdapt) + 1.0, 1.0), -0.5);
        return anneal;
    }

    public static DoubleMatrix defaultTrajectoryLengthBounds(double lowerbound, double upperbound)
    {
        double[] ba = {Math.log(lowerbound), Math.log(upperbound)};
        DoubleMatrix bound = new DoubleMatrix(1, 2, ba);
        return bound;
    }


    public static void main(String args[]) {
//        Random rand = new Random(1);
//        DoubleMatrix targetSigma = new DoubleMatrix(new double[][]
//                {{1.0, 0.99}, {0.99, 1.0}});
//        DoubleMatrix targetMean = new DoubleMatrix( new double[] {3.0, 5.0});
//
//        GaussianExample ge = new ARFS.GaussianExample(targetSigma, targetMean);
//
//        ARFS arfs = initializeARFSWithLBFGS(3000, 1000, ge, ge, 2);
//
//        arfs.keepSamples();
//        arfs.sample(rand);
//        arfs.samples.columnMeans().print();
//
//        arfs.sample(rand);
//        arfs.samples.columnMeans().print();
    }

    public double getFixedTrajectoryLength()
    {
        return fixedTrajectoryLength;
    }


}