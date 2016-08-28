package conifer.rejfreemodels.phylo;

import bayonet.math.EJMLUtils;
import bayonet.math.SparseVector;
import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;
import blang.factors.FactorList;
import blang.variables.RealVariable;
import conifer.ctmc.expfam.CTMCExpFam;
import conifer.ctmc.expfam.CTMCState;
import conifer.ctmc.expfam.ExpFamParameters;
import conifer.local.CollisionContext;
import conifer.local.CollisionFactor;
import conifer.rejfreeutil.StaticUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;
import org.jblas.util.Random;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by crystal on 2016-07-27.
 */
public class TransitionCountFactor implements CollisionFactor {

    @FactorComponent
    public final ExpFamParameters parameters;

    @FactorComponent
    public final FactorList<RealVariable> variables;


    public CTMCExpFam<CTMCState> ctmcExpFam;

    public CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective auxObjective;

    public int state0;
    public int state1;
    public int state1IdxOfBivariateFeatures;


    public TransitionCountFactor(ExpFamParameters parameters, CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective Objective,
                             CTMCExpFam<CTMCState> ctmcExpFam, int state0, int state1,
                                 FactorList<RealVariable>variables, int state1IdxOfBivariateFeatures){
        this.parameters = parameters;
        this.ctmcExpFam = ctmcExpFam;
        this.auxObjective = Objective;
        this.ctmcExpFam = ctmcExpFam;
        this.state0 = state0;
        this.state1 = state1;
        this.variables = variables;
        this.state1IdxOfBivariateFeatures= state1IdxOfBivariateFeatures;

    }



    /**
     * Computer a lower bound for the next collision time.
     *
     * @param context The information relevant to making this calculation
     * @return A pair where the first item is a time, which could be either
     * the collision time (in which case the second item should be true)
     * or a strict lower bound for the collision time (in which case the
     * second item should be false).
     */
    @Override
    public Pair<Double, Boolean> getLowerBoundForCollisionDeltaTime(CollisionContext context) {
        double denominator = 0;
        double result;
        SparseVector [][] bivariateFeatures = ctmcExpFam.bivariateFeatures;
        SparseVector [] univariateFeatures = ctmcExpFam.univariateFeatures;
        final DoubleMatrix v = context.velocity;
        double maxOmega = getOmegaMax(context);
        final double c = StaticUtils.generateUnitRateExponential(context.random);

        if(checkState1InSupportOfState0()){
            denominator = Math.abs(-univariateFeatures[state1].dotProduct(v.toArray())
                    -bivariateFeatures[state0][state1IdxOfBivariateFeatures].dotProduct(v.toArray())+ maxOmega)*getTransitionCount();

        }else{
            denominator = Math.abs(-univariateFeatures[state1].dotProduct(v.toArray())+maxOmega)*getTransitionCount();
        }
        result = c / denominator;
        if(result>0){
            double trueIntensity = getTrueIntensity(context, result);
            double ub = getIntensityUpperBound(context,result);

            //double trueIntensityFromDotProduct = getTrueIntensityUsingDotProduct(context, result);

            double ratio = trueIntensity/ ub;
            Random rand = new Random();
            double V = rand.nextDouble();

            if (V > ratio) {

                return Pair.of(result, false);

            } else {

                return Pair.of(result, true);
            }

        }else{

            return Pair.of(Double.POSITIVE_INFINITY, true);
        }


    }


    public double getTrueIntensityUsingDotProduct(CollisionContext context, double tau){

        double result=0;
        DoubleMatrix v = context.velocity;
        DoubleMatrix gradient = gradientAtNextPosition(context, tau);// since gradient is of logDensity, we should get it with respect to potential energy
        result = gradient.dot(v)*(-1);
        result = Math.max(0, result);
        return result;
    }


    public DoubleMatrix gradientAtNextPosition(CollisionContext context, double tau){

        double [] result = new double [nVariables()];
        SparseVector[][] bivariateFeatures = ctmcExpFam.bivariateFeatures; // S index -> index in support
        SparseVector []  univariateFeatures = ctmcExpFam.univariateFeatures;
        univariateFeatures[state1].linearIncrement( getTransitionCount(),result);

        CTMCExpFam.LearnedReversibleModel ctmcModel = ctmcExpFam.reversibleModelWithParameters(getNextPosition(context,tau));

        for(int state=0; state< ctmcExpFam.nStates; state++){
            univariateFeatures[state].linearIncrement(getTransitionCount()*(-1)*ctmcModel.pi[state], result);
        }
        //check state1 is in support of state0
        if(checkState1InSupportOfState0()){
            bivariateFeatures[state0][state1IdxOfBivariateFeatures].linearIncrement(getTransitionCount(), result);
        }
        return new DoubleMatrix(result);

    }

    public double [] getNextPosition(CollisionContext context, double tau){

        double [] v = context.velocity.toArray();
        double [] delta = new DoubleMatrix(v).mul(tau).toArray();
        DoubleMatrix prevPosition = new DoubleMatrix(getPosition());
        DoubleMatrix deltaDist = new DoubleMatrix(delta);
        DoubleMatrix nextPosition = prevPosition.add(deltaDist);
        return nextPosition.toArray();
    }



    public double getTrueIntensity(CollisionContext context, double tau){
        double result;
        SparseVector [][] bivariateFeatures = ctmcExpFam.bivariateFeatures;
        SparseVector [] univariateFeatures = ctmcExpFam.univariateFeatures;
        double [] v = context.velocity.toArray();

        double part1 = bivariateFeatures[state0][state1IdxOfBivariateFeatures].dotProduct(v);
        double part2 = univariateFeatures[state1].dotProduct(v);

        double [] weightForPi = getPosition();

//        double [] delta = new DoubleMatrix(v).mul(tau).toArray();
//
//        DoubleMatrix prevPosition = new DoubleMatrix(weightForPi);
//        DoubleMatrix deltaDist = new DoubleMatrix(delta);
//        DoubleMatrix nextPosition = prevPosition.add(deltaDist);
//
//        double [] nextPositionArray = nextPosition.toArray();
//
//        CTMCExpFam.LearnedReversibleModel ctmcModel = ctmcExpFam.reversibleModelWithParameters(nextPositionArray);
//        double part3 = 0;
//        double term3;
//        for(int i=0; i< ctmcExpFam.nStates; i++){
//            term3 = ctmcModel.pi[i]*univariateFeatures[i].dotProduct(v);
//
//            part3 = part3+ term3;
//        }
//
//        double result1 = -getTransitionCount()*(part1+part2-part3);

        double numerator = 0;
        double denominator = 0;

        double term1=0;
        double term2=0;

        for(int i = 0; i< ctmcExpFam.nStates; i++){

            term1 = Math.exp(univariateFeatures[i].dotProduct(weightForPi)+univariateFeatures[i].dotProduct(v)*tau)* univariateFeatures[i].dotProduct(v);

            numerator = numerator + term1;
                        }

        for(int i = 0; i < ctmcExpFam.nStates; i++){

            term2 = Math.exp(univariateFeatures[i].dotProduct(weightForPi)+univariateFeatures[i].dotProduct(v)*tau);

            denominator = denominator +  term2;
        }

        result = -getTransitionCount()*(part1+part2-numerator/denominator);

//        if(Math.abs(result-result1)>1e-6)
//            throw new RuntimeException("the two intensities are not equal");

        return Math.max(0, result);
    }


    public double getIntensityUpperBound(CollisionContext context, double tau){
        double result = 0;

        SparseVector [][] bivariateFeatures = ctmcExpFam.bivariateFeatures;
        SparseVector [] univariateFeatures = ctmcExpFam.univariateFeatures;
        double [] v = context.velocity.toArray();
        double part1 = univariateFeatures[state1].dotProduct(v)+bivariateFeatures[state0][state1IdxOfBivariateFeatures].dotProduct(v);
        result = getTransitionCount() * Math.abs(getOmegaMax(context)- part1);
        return result;
    }


    /**
     * @return The gradient of the log-likelihood (NOT the energy, to be
     * in agreement with the inherited Factor)
     */
    @Override
    public DoubleMatrix gradient() {
        double [] result = new double [nVariables()];
        SparseVector[][] bivariateFeatures = ctmcExpFam.bivariateFeatures; // S index -> index in support
        SparseVector []  univariateFeatures = ctmcExpFam.univariateFeatures;
        univariateFeatures[state1].linearIncrement( getTransitionCount(),result);

        CTMCExpFam.LearnedReversibleModel ctmcModel = ctmcExpFam.reversibleModelWithParameters(getPosition());

        for(int state=0; state< ctmcExpFam.nStates; state ++){
            univariateFeatures[state].linearIncrement(getTransitionCount()*(-1)*ctmcModel.pi[state], result);
        }
        //check state1 is in support of state0
        if(checkState1InSupportOfState0()){
            bivariateFeatures[state0][state1IdxOfBivariateFeatures].linearIncrement(getTransitionCount(), result);
        }
        return new DoubleMatrix(result);
    }

    @Override
    public RealVariable getVariable(int gradientCoordinate) {
        return   variables.list.get(gradientCoordinate);
    }

    @Override
    public int nVariables() {
        return parameters.getVector().length;
    }

    /**
     * Note: the density should be normalized, and in log scale.
     *
     * @return The log of the density for the current
     * assignment of parameters and realization.
     */
    @Override
    public double logDensity() {
        return getTransitionCount() * Math.log(getRateMtxElement());
    }

    public double getOmegaMax(CollisionContext context){
        final DoubleMatrix v = context.velocity;
        double maxOmega=0;
        SparseVector[]  univariateFeatures = ctmcExpFam.univariateFeatures;

        double nStates = parameters.globalExponentialFamily.nStates;
        double candidate =0;
        maxOmega = Math.abs(univariateFeatures[0].dotProduct(v.toArray()));
        for(int i=1; i< nStates; i++){
            candidate = Math.abs(univariateFeatures[i].dotProduct(v.toArray()));
            if(candidate > maxOmega)
                maxOmega = candidate;
        }
        return maxOmega;
    }

    public double [] getPosition(){
        double [] result = new double[nVariables()];
        for(int i=0; i< nVariables(); i++){
            result[i] = variables.list.get(i).getValue();
        }
        return result;

    }

    private double getRateMtxElement(){
        CTMCExpFam.LearnedReversibleModel ctmcModel = ctmcExpFam.reversibleModelWithParameters(getPosition());
        return ctmcModel.getRateMatrix()[state0][state1];

    }

    private double getTransitionCount(){

        return auxObjective.nTrans[state0][state1IdxOfBivariateFeatures];
    }

    public boolean checkState1InSupportOfState0(){
        // check if state 1 is in support of state0
        int [][] supports = ctmcExpFam.supports;
        return ArrayUtils.contains(supports[state0], state1);
    }

}
