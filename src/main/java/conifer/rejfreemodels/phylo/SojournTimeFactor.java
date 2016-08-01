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
 * Created by crystal on 2016-07-26.
 */
public class SojournTimeFactor implements CollisionFactor{

    @FactorComponent
    public final ExpFamParameters parameters;

    @FactorArgument(makeStochastic = true)
    public final List<RealVariable> weights;

    public CTMCExpFam<CTMCState> ctmcExpFam;

    public CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective auxObjective;

    public int state0;
    public int state1;

    public SojournTimeFactor(ExpFamParameters parameters, CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective Objective,
                             CTMCExpFam<CTMCState> ctmcExpFam, int state0, int state1, List<RealVariable> weights){
        this.parameters = parameters;
        this.ctmcExpFam = ctmcExpFam;
        this.auxObjective = Objective;
        this.ctmcExpFam = ctmcExpFam;
        this.state0 = state0;
        this.state1 = state1;
        this.weights= weights;
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

        final DoubleMatrix v  = context.velocity;
        double maxOmega = getOmegaMax(context);
        double vdotProdBivariateFeatures = getVelocityDotProductWithBivariateFeature(context);

        double part1 = 1/(2*maxOmega+ vdotProdBivariateFeatures);
        // generate a number where c = -logV, where V is uniformly distributed in (0, 1)
        final double c = StaticUtils.generateUnitRateExponential(context.random);
        double part2 = Math.log(c/(auxObjective.holdTimes[state0]*getRateMtxElement()));
        double result = Double.POSITIVE_INFINITY;
        double candidate = part1 * part2;
        if(candidate>0)
            result = candidate;
        return Pair.of(result, false);

    }

    public double getVelocityDotProductWithBivariateFeature(CollisionContext context){
        final DoubleMatrix v = context.velocity;
        double result;
        SparseVector[][] bivariateFeatures = ctmcExpFam.bivariateFeatures;

        if(checkState1InSupportOfState0()){
            result = bivariateFeatures[state0][state1].dotProduct(v.toArray());
        }else{
            result =0;
        }
        return result;
    }

    public boolean checkState1InSupportOfState0(){
        // check if state 1 is in support of state0
        int [][] supports = ctmcExpFam.supports;
        return ArrayUtils.contains(supports[state0], state1);
    }

    public double getOmegaMax(CollisionContext context){
        final DoubleMatrix v = context.velocity;
        double maxOmega=0;
        SparseVector []  univariateFeatures = ctmcExpFam.univariateFeatures;

        double nStates = ctmcExpFam.nStates;
        double candidate =0;
        maxOmega = Math.abs(univariateFeatures[0].dotProduct(v.toArray()));
        for(int i=1; i< nStates; i++){
            candidate = Math.abs(univariateFeatures[i].dotProduct(v.toArray()));
            if(candidate > maxOmega)
                maxOmega = candidate;
        }
        return maxOmega;
    }

    /**
     * @return The gradient of the log-likelihood (NOT the energy, to be
     * in agreement with the inherited Factor)
     */
    @Override
    public DoubleMatrix gradient() {
        double multiplier= logDensity();
        double [] result = new double [nVariables()];
        SparseVector[][] bivariateFeatures = ctmcExpFam.bivariateFeatures; // S index -> index in support
        SparseVector []  univariateFeatures = ctmcExpFam.univariateFeatures;
        univariateFeatures[state1].linearIncrement( logDensity(),result);

        CTMCExpFam.LearnedReversibleModel ctmcModel = ctmcExpFam.reversibleModelWithParameters(parameters.weights);

        for(int state=0; state< ctmcExpFam.nStates; state ++){
            univariateFeatures[state].linearIncrement(logDensity()*(-1)*ctmcModel.pi[state], result);
        }
        //check state1 is in support of state0
        if(checkState1InSupportOfState0()){
            bivariateFeatures[state0][state1].linearIncrement(logDensity(), result);
        }
        return new DoubleMatrix(result);
    }

    @Override
    public RealVariable getVariable(int gradientCoordinate) {
        return RealVariable.real(parameters.getVector()[gradientCoordinate]);
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
        double result;
        result = getSojournTime()*getRateMtxElement()*(-1.0);
        return result;
    }

    private double getRateMtxElement(){
       return parameters.getRateMatrix(0)[state0][state1];

    }

    private double getSojournTime(){

        return auxObjective.holdTimes[state0];
    }


}
