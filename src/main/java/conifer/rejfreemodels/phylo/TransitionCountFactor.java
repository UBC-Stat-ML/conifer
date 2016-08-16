package conifer.rejfreemodels.phylo;

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

    //@FactorArgument(makeStochastic = true)
    public final List<RealVariable> weights;

    @FactorComponent
    public final FactorList<RealVariable> variables;


    public CTMCExpFam<CTMCState> ctmcExpFam;

    public CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective auxObjective;

    public int state0;
    public int state1;
    public int state1IdxOfBivariateFeatures;


    public TransitionCountFactor(ExpFamParameters parameters, CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective Objective,
                             CTMCExpFam<CTMCState> ctmcExpFam, int state0, int state1,List<RealVariable> weights,
                                 FactorList<RealVariable>variables, int state1IdxOfBivariateFeatures){
        this.parameters = parameters;
        this.ctmcExpFam = ctmcExpFam;
        this.auxObjective = Objective;
        this.ctmcExpFam = ctmcExpFam;
        this.state0 = state0;
        this.state1 = state1;
        this.weights = getWeights();
        this.variables = variables;
        this.state1IdxOfBivariateFeatures= state1IdxOfBivariateFeatures;

    }

    public List<RealVariable> getWeights(){
        return transformWeightIntoReal(parameters);
    }

    public List<RealVariable> transformWeightIntoReal(ExpFamParameters parameters){

        List<RealVariable> result = new ArrayList<>();
        for(int i=0; i<parameters.getVector().length;i++){
            result.add(RealVariable.real(parameters.getVector()[i]));
        }
        return result;
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
            denominator = Math.abs(-univariateFeatures[state1].dotProduct(v.toArray()) +
                    -bivariateFeatures[state0][state1IdxOfBivariateFeatures].dotProduct(v.toArray())+ maxOmega)*getTransitionCount();

        }else{
            denominator = Math.abs(univariateFeatures[state1].dotProduct(v.toArray())+maxOmega)*getTransitionCount();
        }
        result = c / denominator;
        double ratio = getTrueIntensity(context, result)/getIntensityUpperBound(context,result);
        Random rand = new Random();
        double V = rand.nextDouble();

        if(V > ratio){

            return Pair.of(result, false);

        }else{

            return Pair.of(result, true);
        }

    }





    public double getTrueIntensity(CollisionContext context, double tau){
        double result;
        SparseVector [][] bivariateFeatures = ctmcExpFam.bivariateFeatures;
        SparseVector [] univariateFeatures = ctmcExpFam.univariateFeatures;
        double [] v = context.velocity.toArray();

        double part1 = bivariateFeatures[state0][state1IdxOfBivariateFeatures].dotProduct(v);
        double part2 = univariateFeatures[state1].dotProduct(v);

        double [] weightForPi;
        CTMCExpFam.LearnedReversibleModel ctmcModel = ctmcExpFam.reversibleModelWithParameters(parameters.weights);
        weightForPi = ctmcModel.weights;

        double numerator = 0;
        double denominator = 0;

        for(int i = 0; i< ctmcExpFam.nStates; i++){

            numerator = numerator +
                    Math.exp(univariateFeatures[i].dotProduct(weightForPi))* Math.exp(univariateFeatures[i].dotProduct(v)*tau)* univariateFeatures[i].dotProduct(v);
        }

        for(int i = 0; i < ctmcExpFam.nStates; i++){

            denominator = denominator + Math.exp(univariateFeatures[i].dotProduct(weightForPi))* Math.exp(univariateFeatures[i].dotProduct(v)*tau);
        }

        result = -getTransitionCount()*(part1+part2-numerator/denominator);
        return Math.max(0, result);
    }


    public double getIntensityUpperBound(CollisionContext context, double tau){
        double result = 0;

        SparseVector [][] bivariateFeatures = ctmcExpFam.bivariateFeatures;
        SparseVector [] univariateFeatures = ctmcExpFam.univariateFeatures;
        double [] v = context.velocity.toArray();
        double part1 = univariateFeatures[state1].dotProduct(v)+bivariateFeatures[state0][state1IdxOfBivariateFeatures].dotProduct(v);
        result = getTransitionCount() * Math.abs(getOmegaMax(context)- part1);
        result = Math.max(0, result);
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

        CTMCExpFam.LearnedReversibleModel ctmcModel = ctmcExpFam.reversibleModelWithParameters(parameters.weights);

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

    private double getRateMtxElement(){
        return parameters.getRateMatrix(0)[state0][state1];

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
