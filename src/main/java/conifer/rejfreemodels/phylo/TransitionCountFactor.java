package conifer.rejfreemodels.phylo;

import bayonet.math.SparseVector;
import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;
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
import java.util.ArrayList;
import java.util.List;

/**
 * Created by crystal on 2016-07-27.
 */
public class TransitionCountFactor implements CollisionFactor {

    @FactorComponent
    public final ExpFamParameters parameters;

   @FactorArgument(makeStochastic = true)
   public final List<RealVariable> weights;


    public CTMCExpFam<CTMCState> ctmcExpFam;

    public CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective auxObjective;

    public int state0;
    public int state1;


    public TransitionCountFactor(ExpFamParameters parameters, CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective Objective,
                             CTMCExpFam<CTMCState> ctmcExpFam, int state0, int state1,List<RealVariable> weights){
        this.parameters = parameters;
        this.ctmcExpFam = ctmcExpFam;
        this.auxObjective = Objective;
        this.ctmcExpFam = ctmcExpFam;
        this.state0 = state0;
        this.state1 = state1;
        this.weights = getWeights();

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
            denominator = Math.abs(univariateFeatures[state1].dotProduct(v.toArray()) +
                    bivariateFeatures[state0][state1].dotProduct(v.toArray())+ maxOmega)*getTransitionCount();

        }else{
            denominator = Math.abs(univariateFeatures[state1].dotProduct(v.toArray())+maxOmega)*getTransitionCount();
        }
        result = c / denominator;
        return Pair.of(result, false);
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
            bivariateFeatures[state0][state1].linearIncrement(getTransitionCount(), result);
        }
        return new DoubleMatrix(result);
    }

    @Override
    public RealVariable getVariable(int gradientCoordinate) {
        return  RealVariable.real(parameters.getVector()[gradientCoordinate]);
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

        return auxObjective.nTrans[state0][state1];
    }

    public boolean checkState1InSupportOfState0(){
        // check if state 1 is in support of state0
        int [][] supports = ctmcExpFam.supports;
        return ArrayUtils.contains(supports[state0], state1);
    }

}
