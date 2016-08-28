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
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.StatUtils;
import org.jblas.DoubleMatrix;
import org.jblas.util.Random;

/**
 * Created by crystal on 2016-07-26.
 */
public class SojournTimeFactor implements CollisionFactor{

    @FactorComponent
    public final ExpFamParameters parameters;


    public CTMCExpFam<CTMCState> ctmcExpFam;

    public CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective auxObjective;

    @FactorComponent
    public final FactorList<RealVariable> variables;

    public int state0;
    public int state1;
    public int state1IdxOfBivariateFeatures;

    public SojournTimeFactor(ExpFamParameters parameters, CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective Objective,
                             CTMCExpFam<CTMCState> ctmcExpFam, int state0, int state1,
                             FactorList<RealVariable> variables, int state1IdxOfBivariateFeatures){
        this.parameters = parameters;
        this.ctmcExpFam = ctmcExpFam;
        this.auxObjective = Objective;
        this.ctmcExpFam = ctmcExpFam;
        this.state0 = state0;
        this.state1 = state1;
        this.variables = variables;
        this.state1IdxOfBivariateFeatures=state1IdxOfBivariateFeatures;

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

        final double [] v  = context.velocity.toArray();
        double maxOmega = getOmegaMax(context);
        double vdotProdBivariateFeatures = getVelocityDotProductWithBivariateFeature(context);

        SparseVector [] univariateFeatures = ctmcExpFam.univariateFeatures;
        SparseVector [][] bivariateFeatures = ctmcExpFam.bivariateFeatures;

        double c = StaticUtils.generateUnitRateExponential(context.random);

        double condition = univariateFeatures[state1].dotProduct(v)+ vdotProdBivariateFeatures;
        double denom = maxOmega + condition;
        double result = 0;

        if(condition >= 0){

            result = 1/(denom) * Math.log(1 + c/(getRateMtxElement()*getSojournTime()));

        }else{

            double term2 = denom/(Math.abs(condition) + maxOmega);
            result = 1/denom * Math.log(1 + c/(getRateMtxElement()*getSojournTime())* term2 );
        }

        double trueIntensity = getTrueIntensity(context, result);
        double ub = getIntensityUpperBound(context,result);

        // once ensure the correctness, we need to comment this
//        double trueIntensityFromDotProduct = getTrueIntensityUsingDotProduct(context, result);
//        if(Math.abs(trueIntensity-trueIntensityFromDotProduct)>1e-6){
//            System.out.print(Math.abs(trueIntensity-trueIntensityFromDotProduct));
//            throw new RuntimeException("the two intensities obtained from two different methods are different");
//        }


        double ratio = trueIntensity/ ub;
        Random rand = new Random();
        double V = rand.nextDouble();
        if(V > ratio){

            return Pair.of(result, false);

        }else{

            return Pair.of(result, true);
        }

    }

    public double [] getPosition(){
        double [] result = new double[nVariables()];
        for(int i=0; i< nVariables(); i++){
            result[i] = variables.list.get(i).getValue();
        }
        return result;

    }

    public double getVelocityDotProductWithBivariateFeature(CollisionContext context){
        final DoubleMatrix v = context.velocity;
        double result;
        SparseVector[][] bivariateFeatures = ctmcExpFam.bivariateFeatures;

        if(checkState1InSupportOfState0()){
            result = bivariateFeatures[state0][state1IdxOfBivariateFeatures].dotProduct(v.toArray());
        }else{
            result =0;
        }
        return result;
    }

    public double getUnivariateFeatureNormalization(){

        SparseVector []  univariateFeatures = ctmcExpFam.univariateFeatures;
        int nStates = ctmcExpFam.nStates;
        double [] pi = new double[nStates];

        double [] weightForPi = getPosition();
        CTMCExpFam.LearnedReversibleModel ctmcModel = ctmcExpFam.reversibleModelWithParameters(weightForPi);
        for (int i = 0; i < nStates; i++)
            pi[i] = Math.exp(univariateFeatures[i].dotProduct(weightForPi));

        double normalization;
        // get the sum of all elements in pi
        normalization = StatUtils.sum(pi);
        return normalization;
    }

    public double [] getNextPosition(CollisionContext context, double tau){

        double [] v = context.velocity.toArray();
        double [] delta = new DoubleMatrix(v).mul(tau).toArray();
        DoubleMatrix prevPosition = new DoubleMatrix(getPosition());
        DoubleMatrix deltaDist = new DoubleMatrix(delta);
        DoubleMatrix nextPosition = prevPosition.add(deltaDist);
        return nextPosition.toArray();
    }

    public DoubleMatrix gradientAtNextPosition(CollisionContext context, double tau) {

        double [] result = new double [nVariables()];
        SparseVector[][] bivariateFeatures = ctmcExpFam.bivariateFeatures; // S index -> index in support
        SparseVector []  univariateFeatures = ctmcExpFam.univariateFeatures;
        CTMCExpFam.LearnedReversibleModel ctmcModel = ctmcExpFam.reversibleModelWithParameters(getNextPosition(context, tau));
        double rateMtxElement = ctmcModel.getRateMatrix()[state0][state1];
        double newLogDensity = getSojournTime()*rateMtxElement*(-1.0);

        univariateFeatures[state1].linearIncrement( newLogDensity,result);

        for(int state=0; state< ctmcExpFam.nStates; state ++){
            univariateFeatures[state].linearIncrement(newLogDensity*(-1)*ctmcModel.pi[state], result);
        }
        //check state1 is in support of state0
        if(checkState1InSupportOfState0()){
            bivariateFeatures[state0][state1IdxOfBivariateFeatures].linearIncrement(newLogDensity, result);
        }
        return new DoubleMatrix(result);
    }

    public double getTrueIntensityUsingDotProduct(CollisionContext context, double tau){

        double result=0;
        DoubleMatrix v = context.velocity;
        DoubleMatrix gradient = gradientAtNextPosition(context, tau);// since gradient is of logDensity, we should get it with respect to potential energy
        result = gradient.dot(v)*(-1);
        result = Math.max(0, result);
        return result;
    }

    public double getTrueIntensity(CollisionContext context, double tau){

        double piUniformization = getUnivariateFeatureNormalization();
        double sojournTime = getSojournTime();
        double rateMtxEntry = getRateMtxElement();
        double part1 = piUniformization*sojournTime*rateMtxEntry;
        SparseVector []  univariateFeatures = ctmcExpFam.univariateFeatures;
        SparseVector[][] bivariateFeatures = ctmcExpFam.bivariateFeatures;
        double [] v= context.velocity.toArray();
        double part2 = (univariateFeatures[state1].dotProduct(v)+bivariateFeatures[state0][state1IdxOfBivariateFeatures].dotProduct(v));
        double part3 = Math.exp(part2*tau);

        double [] weightForPi = getPosition();

        double denominator1 = 0;
        for(int i=0; i<ctmcExpFam.nStates;i++){
            denominator1 = denominator1 + Math.exp(univariateFeatures[i].dotProduct(getNextPosition(context, tau)));
        }

        double part4 = part2/denominator1;

        double denominator2 = denominator1*denominator1;
        double numerator2 =0;
        for(int i=0; i<ctmcExpFam.nStates;i++){
            numerator2 = numerator2 + Math.exp(univariateFeatures[i].dotProduct(getNextPosition(context,tau)))*univariateFeatures[i].dotProduct(v);
        }

        double part5 = numerator2/denominator2;
        double result = 0;
        result = part1 * part3 * (part4 -part5);
        result = Math.max(0, result);

        return result;
    }


    public double getIntensityUpperBound(CollisionContext context, double tau){
        double result = 0;
        SparseVector[][] bivariateFeatures = ctmcExpFam.bivariateFeatures;
        SparseVector []  univariateFeatures = ctmcExpFam.univariateFeatures;
        double [] v = context.velocity.toArray();
        double maxOmega = getOmegaMax(context);

        double condition = univariateFeatures[state1].dotProduct(v)+bivariateFeatures[state0][state1IdxOfBivariateFeatures].dotProduct(v);
        double term1 = 0;
        if(condition >= 0){
            term1 = (condition+maxOmega);
            result = getSojournTime()* getRateMtxElement() * Math.exp(term1 * tau) * term1;
        }else{
            term1 = Math.abs(condition)+ maxOmega;
            result = getSojournTime()* getRateMtxElement() * Math.exp((condition + maxOmega)* tau) * term1;
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
        double [] result = new double [nVariables()];
        SparseVector[][] bivariateFeatures = ctmcExpFam.bivariateFeatures; // S index -> index in support
        SparseVector []  univariateFeatures = ctmcExpFam.univariateFeatures;
        univariateFeatures[state1].linearIncrement( logDensity(),result);

        CTMCExpFam.LearnedReversibleModel ctmcModel = ctmcExpFam.reversibleModelWithParameters(getPosition());

        for(int state=0; state< ctmcExpFam.nStates; state ++){
            univariateFeatures[state].linearIncrement(logDensity()*(-1)*ctmcModel.pi[state], result);
        }
        //check state1 is in support of state0
        if(checkState1InSupportOfState0()){
            bivariateFeatures[state0][state1IdxOfBivariateFeatures].linearIncrement(logDensity(), result);
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
        CTMCExpFam.LearnedReversibleModel ctmcModel = ctmcExpFam.reversibleModelWithParameters(getPosition());
        return ctmcModel.getRateMatrix()[state0][state1];

    }

    private double getSojournTime(){

        return auxObjective.holdTimes[state0];
    }


}
