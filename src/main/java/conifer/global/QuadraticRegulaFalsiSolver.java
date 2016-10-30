package conifer.global;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.*;
import org.apache.commons.math3.exception.*;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Incrementor;
import org.apache.commons.math3.util.MathUtils;
import org.jblas.DoubleMatrix;

import conifer.global.GlobalRFSampler.CollisionSolver;
import bayonet.math.NumericalUtils;
import bayonet.opt.DifferentiableFunction;
import bayonet.opt.LBFGSMinimizer;


/**
 * Created by crystal on 2016-10-27.
 */
public class QuadraticRegulaFalsiSolver {

    /** Default relative accuracy. */
    private static final double DEFAULT_RELATIVE_ACCURACY = 1e-14;
    /** Default function value accuracy. */
    private static final double DEFAULT_FUNCTION_VALUE_ACCURACY = 1e-15;

    protected static final double DEFAULT_ABSOLUTE_ACCURACY = 1e-6;
    private double absoluteAccuracy = DEFAULT_ABSOLUTE_ACCURACY;
    private AllowedSolution allowed;

    private double searchMin;
    /** Higher end of search interval. */
    private double searchMax;
    /** Initial guess. */
    private double searchStart;
    /** Function to solve. */
    private UnivariateFunction function;

    private double  relativeAccuracy = DEFAULT_RELATIVE_ACCURACY;
    private double  functionValueAccuracy = DEFAULT_FUNCTION_VALUE_ACCURACY;


    /** Evaluations counter. */
    private final Incrementor evaluations = new Incrementor();


    protected QuadraticRegulaFalsiSolver(){
        this.absoluteAccuracy = DEFAULT_ABSOLUTE_ACCURACY;
        this.allowed = AllowedSolution.ABOVE_SIDE;

    }

    /**
     * Construct a solver.
     *
     * @param absoluteAccuracy Absolute accuracy.
     */
    protected QuadraticRegulaFalsiSolver(final double absoluteAccuracy) {
        this.absoluteAccuracy = absoluteAccuracy;
        this.allowed = AllowedSolution.ANY_SIDE;
    }

    protected QuadraticRegulaFalsiSolver(final double relativeAccuracy,
                               final double absoluteAccuracy) {
        this.relativeAccuracy = relativeAccuracy;
        this.absoluteAccuracy = absoluteAccuracy;
        this.allowed = AllowedSolution.ANY_SIDE;
    }

    /**
     * Construct a solver.
     *
     * @param relativeAccuracy Maximum relative error.
     * @param absoluteAccuracy Maximum absolute error.
     * @param functionValueAccuracy Maximum function value error.
     */
    protected QuadraticRegulaFalsiSolver(final double relativeAccuracy, final double absoluteAccuracy,
                               final double functionValueAccuracy) {

        this.relativeAccuracy = relativeAccuracy;
        this.absoluteAccuracy = absoluteAccuracy;
        this.functionValueAccuracy = functionValueAccuracy;
        this.allowed = AllowedSolution.ANY_SIDE;
    }

    public double getAbsoluteAccuracy() {
        return absoluteAccuracy;
    }
    /**
     * {@inheritDoc}
     */
    public double getRelativeAccuracy() {
        return relativeAccuracy;
    }
    /**
     * {@inheritDoc}
     */
    public double getFunctionValueAccuracy() {
        return functionValueAccuracy;
    }



    /**
     * @return the lower end of the search interval.
     */
    public double getMin() {
        return searchMin;
    }
    /**
     * @return the higher end of the search interval.
     */
    public double getMax() {

        return searchMax;
    }

    public double solve(final int maxEval, final UnivariateFunction f,
                        final double min, final double max,
                        final AllowedSolution allowedSolution) {
        return solve(maxEval, f, min, max, min + 0.5 * (max - min), allowedSolution);
    }

    public double solve(final int maxEval, final UnivariateFunction f, final double min, final double max){

        return solve(maxEval, f, min, max, min + 0.5 * (max - min), AllowedSolution.ANY_SIDE);
    }

    public double solve(final int maxEval, final UnivariateFunction f,
                        final double min, final double max, final double startValue,
                        final AllowedSolution allowedSolution) {
        this.allowed = allowedSolution;
        return solve(maxEval, f, min, max, startValue);
    }



    public double solve(int maxEval, UnivariateFunction f, double min, double max, double startValue)
            throws TooManyEvaluationsException,
            NoBracketingException {
        // Initialization.
        setup(maxEval, f, min, max, startValue);

        // Perform computation.
        return doSolve();
    }

    protected double computeObjectiveValue(double point)
            throws TooManyEvaluationsException {
        //incrementEvaluationCount();
        return function.value(point);
    }


    protected final double getHn(double x0, double x1, double f0, double f1){
        double hn;
        hn = (x1-x0)/(f1-f0);
        return hn;
    }


    protected final double getCn(double x0, double x1, double f0, double f1){
        double cn;
        cn = (f1*x0-f0*x1)/(f1-f0);
        return cn;
    }


    protected final double getPn(double fxn, double fcn){
        double pn = Math.signum(fxn - fcn);
        return pn;
    }

    protected final double doSolve()
            throws ConvergenceException {
        // Get initial solution
        double x0 = getMin();
        double x1 = getMax();
        double f0 = computeObjectiveValue(x0);
        double f1 = computeObjectiveValue(x1);

        // If one of the bounds is the exact root, return it. Since these are
        // not under-approximations or over-approximations, we can return them
        // regardless of the allowed solutions.
        if (f0 == 0.0) {
            return x0;
        }
        if (f1 == 0.0) {
            return x1;
        }

        // Verify bracketing of initial solution.
        verifyBracketing(x0, x1);

        // Get accuracies.
        final double ftol = getFunctionValueAccuracy();
        final double atol = getAbsoluteAccuracy();
        final double rtol = getRelativeAccuracy();
        double anBar = x0;
        double bnBar = x1;

        double xn = searchStart;
        double xStar = xn;
        double pn, hn, cn = 0;
        // Keep track of inverted intervals, meaning that the left bound is
        // larger than the right bound.
        boolean inverted = false;

        double fanBar = 0;
        double fbnBar = 0;
        double fxn = computeObjectiveValue(xn);

        // Keep finding better approximations.
        while (true) {
            //get hn, pn, cn first and then get cnbar
            hn = getHn(x0, x1, f0, f1);
            cn = getCn(x0, x1, f0, f1);
            final double fcn = computeObjectiveValue(cn);
            pn = getPn(fxn, fcn);

            if(fcn == 0|| FastMath.abs(fcn)<= ftol){
                return cn;
            }

            // Update the bounds with the new approximation.
            if (f0 * fcn < 0) {
                // The value of x1 has switched to the other bound, thus inverting
                // the interval.
                anBar = x0;
                fanBar = f0;
                bnBar = cn;
                fbnBar = fcn;

            }

            if (f1 * fcn < 0){
                anBar = cn;
                fanBar = fcn;
                bnBar = x1;
                fbnBar = f1;

            }

            // check convergence of cn sequence

            if (FastMath.abs(cn - xStar) < FastMath.max(rtol * FastMath.abs(cn), atol)) {
                return cn;
            }

            xStar = cn;

            // call exponential iterative procedure
            // Calculate the next approximation.
            final double cnBar = xn * Math.exp(-hn*fxn*fxn/(xn*(pn*fxn*fxn+fxn-fcn)));
            final double fcnBar = computeObjectiveValue(cnBar);

            incrementEvaluationCount();



            // If the new approximation is the exact root, return it. Since
            // this is not an under-approximation or an over-approximation,
            // we can return it regardless of the allowed solutions.
            if (fcnBar == 0.0) {
                return cnBar;
            }


            if(fanBar * fcnBar < 0){
                x0 = anBar;
                x1 = cnBar;
                f0 = fanBar;
                f1 = fcnBar;


            }else{
                x0 = cnBar;
                f0 = fcnBar;
                x1 = bnBar;
                f1 = fbnBar;
            }

            if(FastMath.abs(cnBar-xn) < FastMath.max(rtol * FastMath.abs(Math.max(cnBar, xn)),
                    atol)){
                return xn;
            }

            xn = cnBar;
            fxn = fcnBar;

            // If the function value of the last approximation is too small,
            // given the function value accuracy, then we can't get closer to
            // the root than we already are.
            if (FastMath.abs(fxn) <= ftol) {

                return xn;
            }

            // If the current interval is within the given accuracies, we
            // are satisfied with the current approximation.
            if (FastMath.abs(x1 - x0) < FastMath.max(rtol * FastMath.abs(x1),
                    atol)) {
                switch (allowed) {
                    case ANY_SIDE:
                        return x1;
                    case LEFT_SIDE:
                        return inverted ? x1 : x0;
                    case RIGHT_SIDE:
                        return inverted ? x0 : x1;
                    case BELOW_SIDE:
                        return (f1 <= 0) ? x1 : x0;
                    case ABOVE_SIDE:
                        return (f1 >= 0) ? x1 : x0;
                    default:
                        throw new MathInternalError();
                }
            }



            if(cnBar < anBar || cnBar > bnBar){
                x0 = anBar;
                f0 = fanBar;
                x1 = bnBar;
                f1 = fbnBar;
                if(cnBar < anBar){
                    xn = anBar;
                    fxn = fanBar;
                }
                if(cnBar >bnBar){
                    xn = bnBar;
                    fxn = fbnBar;
                }
            }

        }
    }

    /**
     * Check that the endpoints specify an interval and the function takes
     * opposite signs at the endpoints.
     *
     * @param lower Lower endpoint.
     * @param upper Upper endpoint.
     * @throws NullArgumentException if the function has not been set.
     * @throws NoBracketingException if the function has the same sign at
     * the endpoints.
     */
    protected void verifyBracketing(final double lower,
                                    final double upper)
            throws NullArgumentException,
            NoBracketingException {
        UnivariateSolverUtils.verifyBracketing(function, lower, upper);
    }


    /**
     * Prepare for computation.
     * Subclasses must call this method if they override any of the
     * {@code solve} methods.
     *
     * @param f Function to solve.
     * @param min Lower bound for the interval.
     * @param max Upper bound for the interval.
     * @param startValue Start value to use.
     * @param maxEval Maximum number of evaluations.
     * @exception NullArgumentException if f is null
     */
    protected void setup(int maxEval,
                         UnivariateFunction f,
                         double min, double max,
                         double startValue)
            throws NullArgumentException {
        // Checks.
        MathUtils.checkNotNull(f);

        // Reset.
        searchMin = min;
        searchMax = max;
        searchStart = startValue;
        function = f;
        evaluations.setMaximalCount(maxEval);
        evaluations.resetCount();
    }

    /**
     * Increment the evaluation count by one.
     * Method {@link #computeObjectiveValue(double)} calls this method internally.
     * It is provided for subclasses that do not exclusively use
     * {@code computeObjectiveValue} to solve the function.
     * See e.g. {@link AbstractUnivariateDifferentiableSolver}.
     *
     * @throws TooManyEvaluationsException when the allowed number of function
     * evaluations has been exhausted.
     */
    protected void incrementEvaluationCount()
            throws TooManyEvaluationsException {
        try {
            evaluations.incrementCount();
        } catch (MaxCountExceededException e) {
            throw new TooManyEvaluationsException(e.getMax());
        }
    }
}
