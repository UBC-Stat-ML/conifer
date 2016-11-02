package conifer.global;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.AbstractUnivariateDifferentiableSolver;
import org.apache.commons.math3.analysis.solvers.AllowedSolution;
import org.apache.commons.math3.analysis.solvers.UnivariateSolverUtils;
import org.apache.commons.math3.exception.*;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Incrementor;
import org.apache.commons.math3.util.MathUtils;

/**
 * Created by crystal on 2016-11-01.
 */
public abstract class SelfImplementedSolver {
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


    protected SelfImplementedSolver(){
        this.absoluteAccuracy = DEFAULT_ABSOLUTE_ACCURACY;
        this.allowed = AllowedSolution.ABOVE_SIDE;

    }

    /**
     * Construct a solver.
     *
     * @param absoluteAccuracy Absolute accuracy.
     */
    protected SelfImplementedSolver(final double absoluteAccuracy) {
        this.absoluteAccuracy = absoluteAccuracy;
        this.allowed = AllowedSolution.ANY_SIDE;
    }

    protected SelfImplementedSolver(final double relativeAccuracy,
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
    protected SelfImplementedSolver(final double relativeAccuracy, final double absoluteAccuracy,
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



    protected final double getCn(double x0, double x1, double f0, double f1){
        double cn;
        cn = (f1*x0-f0*x1)/(f1-f0);
        return cn;
    }

    protected abstract double doSolve()
            throws TooManyEvaluationsException, NoBracketingException;



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
