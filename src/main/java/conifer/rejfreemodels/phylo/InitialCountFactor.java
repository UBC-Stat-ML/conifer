package conifer.rejfreemodels.phylo;

import bayonet.math.SparseVector;
import blang.annotations.FactorComponent;
import blang.factors.FactorList;
import blang.variables.RealVariable;
import conifer.ctmc.expfam.CTMCExpFam;
import conifer.ctmc.expfam.CTMCState;
import conifer.ctmc.expfam.ExpFamParameters;
import conifer.local.CollisionContext;
import conifer.local.CollisionFactor;
import conifer.rejfreeutil.StaticUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;
import org.jblas.util.Random;

/**
 * Created by crystal on 2016-07-28.
 */
public class InitialCountFactor implements CollisionFactor {

    @FactorComponent
    public final ExpFamParameters parameters;


    @FactorComponent
    public final FactorList<RealVariable> variables;

    public CTMCExpFam<CTMCState> ctmcExpFam;

    public CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective auxObjective;

    public int state0;


    public InitialCountFactor(ExpFamParameters parameters, CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective Objective,
                                 CTMCExpFam<CTMCState> ctmcExpFam, int state0, FactorList<RealVariable> variables){
        this.parameters = parameters;
        this.ctmcExpFam = ctmcExpFam;
        this.auxObjective = Objective;
        this.ctmcExpFam = ctmcExpFam;
        this.state0 = state0;
        this.variables = variables;
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
        double result;
        double denominator;
        SparseVector [] univariateFeatures = ctmcExpFam.univariateFeatures;
        final DoubleMatrix v = context.velocity;
        double maxOmega = getOmegaMax(context);
        final double c = StaticUtils.generateUnitRateExponential(context.random);
        denominator = maxOmega-univariateFeatures[state0].dotProduct(v.toArray());

        //check if denominator is very close to zero; if so, result is positive infinity
        if (denominator < 2 * Double.MIN_VALUE){
            result = Double.POSITIVE_INFINITY;
            return Pair.of(result, false);
        }else{
            result = c/(denominator*getInitialCount());
            double trueIntensity = getTrueIntensity(context, result);
            double ub = getIntensityUpperBound(context, result);
            double ratio = trueIntensity/ub;

            Random rand = new Random();
            double V = rand.nextDouble();
            if(V > ratio){
                return Pair.of(result, false);
            }else{
                return Pair.of(result, true);
            }
        }


    }

    public double [] getPosition(){
        double [] result = new double[nVariables()];
        for(int i=0; i< nVariables(); i++){
            result[i] = variables.list.get(i).getValue();
        }
        return result;

    }

    public double getTrueIntensity(CollisionContext context, double tau){
        SparseVector [] univariateFeatures = ctmcExpFam.univariateFeatures;
        double []  v = context.velocity.toArray();

        double part1 = univariateFeatures[state0].dotProduct(v);
        double numerator = 0;
        double denominator=0;

        double [] weightForPi = getPosition();


        for(int i=0; i < ctmcExpFam.nStates; i++){
            numerator = numerator + Math.exp(univariateFeatures[i].dotProduct(v)*tau + univariateFeatures[i].dotProduct(weightForPi))*univariateFeatures[i].dotProduct(v);
        }

        for(int i=0; i < ctmcExpFam.nStates; i++){
            denominator = denominator + Math.exp(univariateFeatures[i].dotProduct(v)*tau + univariateFeatures[i].dotProduct(weightForPi));
        }

        double result;
        result = -getInitialCount()*(part1-numerator/denominator);
        result = Math.max(0, result);
        return result;
    }



    public double getIntensityUpperBound(CollisionContext context, double tau){

        SparseVector [] univariateFeatures = ctmcExpFam.univariateFeatures;
        final DoubleMatrix v = context.velocity;
        double maxOmega = getOmegaMax(context);
        double result;
        result = getInitialCount()*(maxOmega-univariateFeatures[state0].dotProduct(v.toArray()));
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
        SparseVector []  univariateFeatures = ctmcExpFam.univariateFeatures;
        univariateFeatures[state0].linearIncrement( getInitialCount(),result);

        CTMCExpFam.LearnedReversibleModel ctmcModel = ctmcExpFam.reversibleModelWithParameters(getPosition());

        for(int state=0; state< ctmcExpFam.nStates; state ++){
            univariateFeatures[state].linearIncrement(getInitialCount()*(-1)*ctmcModel.pi[state], result);
        }
        return new DoubleMatrix(result);
    }

    @Override
    public RealVariable getVariable(int gradientCoordinate) {
        return  variables.list.get(gradientCoordinate);
    }

    @Override
    public int nVariables() {
        return parameters.getVector().length;
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

    /**
     * Note: the density should be normalized, and in log scale.
     *
     * @return The log of the density for the current
     * assignment of parameters and realization.
     */
    @Override
    public double logDensity() {

        CTMCExpFam.LearnedReversibleModel ctmcModel = ctmcExpFam.reversibleModelWithParameters(getPosition());
        return getInitialCount() * Math.log(ctmcModel.pi[state0]);
    }

    public double getInitialCount(){

        return auxObjective.nInit[state0];

    }

}
