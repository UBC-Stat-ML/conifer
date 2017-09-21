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
 * This solver implements the algorithm of the paper
 * "AN IMPROVED PEGASUS METHOD FOR ROOT FINDING " by Richard F. King
 * Created by crystal on 2016-10-31.
 */
public class ImprovedPegasusSolver extends SelfImplementedSolver{

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


    protected ImprovedPegasusSolver(){
        super();

    }

    /**
     * Construct a solver.
     *
     * @param absoluteAccuracy Absolute accuracy.
     */
    protected ImprovedPegasusSolver(final double absoluteAccuracy) {
        super(absoluteAccuracy);
    }

    protected ImprovedPegasusSolver(final double relativeAccuracy,
                                         final double absoluteAccuracy) {
        super(relativeAccuracy, absoluteAccuracy);
    }

    /**
     * Construct a solver.
     *
     * @param relativeAccuracy Maximum relative error.
     * @param absoluteAccuracy Maximum absolute error.
     * @param functionValueAccuracy Maximum function value error.
     */
    protected ImprovedPegasusSolver(final double relativeAccuracy, final double absoluteAccuracy,
                                         final double functionValueAccuracy) {

        super(relativeAccuracy, absoluteAccuracy, functionValueAccuracy);

    }


    @Override
    protected final double doSolve()
            throws ConvergenceException {
        // Get initial solution
        double x0 = getMin();
        double x1 = getMax();
        double f0 = computeObjectiveValue(x0);
        double f1 = computeObjectiveValue(x1);
        double tmp = 0;
        double ftmp = 0;

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

        double xn = searchStart;
        double xStar = xn;

        // Keep finding better approximations.
        while (true) {

            // step 3
            xn = getCn(x0, x1, f0, f1);
            double fxn = computeObjectiveValue(xn);
            incrementEvaluationCount();

            if(fxn == 0|| FastMath.abs(fxn)<= ftol){
                return xn;
            }

            // check convergence of xn sequence
            if (FastMath.abs(xn - xStar) < FastMath.max(rtol * FastMath.abs(xn), atol)) {
                return xn;
            }

            // step 4 according to the paper interchange (x0, f0) and (x1, f1) if f1fxn <0
            if(f1 * fxn < 0){
                tmp = x0;
                ftmp = f0;
                x0 = x1;
                f0 = f1;
                x1 = tmp;
                f1 = ftmp;
            }

            // step 5, with the latest (x1, f1) (x2, f2) if (f1f2 > 0) replace (x0, f0)
            // by (x0, f0f1/(f1+f2)) and (x1, f1) by (x2, f2) use step 3 to get new (x2, f2)
            if(f1 * fxn > 0){
                f0 = f0 * f1/(f1 + fxn);
                x1 = xn;
                f1 = fxn;
                xn = getCn(x0, x1, f0, f1);
                fxn = computeObjectiveValue(xn);
                incrementEvaluationCount();
            }

            // step 6 if f1 *fxn < 0, replace (x0, f0) by (x1, f1) and  (x1, f1) by (x2, f2)
            // and go to step 3. Otherwise, go to step 5.
            if(f1 * fxn < 0){
                x0 = x1;
                f0 = f1;
                x1 = xn;
                f1 = fxn;
            }else{
                while(f1 * fxn > 0){
                    f0 = f0* f1/(f1 + fxn);
                    x1 = xn;
                    f1 = fxn;
                    xn = getCn(x0, x1, f0, f1);
                    fxn = computeObjectiveValue(xn);
                    incrementEvaluationCount();
                }
                if(f1 * fxn < 0) {
                    x0 = x1;
                    f0 = f1;
                    x1 = xn;
                    f1 = fxn;
                }

            }

            xStar = xn;

            // If the function value of the last approximation is too small,
            // given the function value accuracy, then we can't get closer to
            // the root than we already are.
            if (fxn ==0 || FastMath.abs(fxn) <= ftol) {

                return xn;
            }

            // If the current interval is within the given accuracies, we
            // are satisfied with the current approximation.
            if (FastMath.abs(x1 - x0) < FastMath.max(rtol * FastMath.abs(x1),
                    atol)) {
                return x1;
            }

        }
    }

}
