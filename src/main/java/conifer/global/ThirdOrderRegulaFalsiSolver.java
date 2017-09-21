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
 * Created by crystal on 2016-10-29.
 */
public class ThirdOrderRegulaFalsiSolver extends SelfImplementedSolver {
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


    protected ThirdOrderRegulaFalsiSolver(){
        super();

    }

    /**
     * Construct a solver.
     *
     * @param absoluteAccuracy Absolute accuracy.
     */
    protected ThirdOrderRegulaFalsiSolver(final double absoluteAccuracy) {
       super(absoluteAccuracy);
    }

    protected ThirdOrderRegulaFalsiSolver(final double relativeAccuracy,
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
    protected ThirdOrderRegulaFalsiSolver(final double relativeAccuracy, final double absoluteAccuracy,
                                         final double functionValueAccuracy) {
        super(relativeAccuracy, absoluteAccuracy, functionValueAccuracy);
    }



    protected final double getQn(double fxn, double f1, double f0){
        double qn = FastMath.abs(fxn)/(f1 - f0);
        return qn;
    }

    protected final double getHn(double x1, double x0, double fxn, double qn){

        double hn = (x1-x0)*qn/(FastMath.abs(fxn));
        return hn;
    }



    protected final double getPn(double xn, double fxn, double hn){
        double input1 = xn - fxn;
        double finput1 = computeObjectiveValue(input1);
        double input2 = xn + fxn;
        double finput2 = computeObjectiveValue(input2);
        double input3 = xn - hn * fxn;
        double finput3 = computeObjectiveValue(input3);

        double numerator = finput3 * (finput1 + finput2 -2 * fxn);
        double denominator = 2*(fxn - finput1)* fxn * fxn;

        double pn = -hn * (numerator / denominator + 1/ (2* xn));
        return pn;
    }

    @Override
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
        double fanBar = f0;
        double fbnBar = f1;

        //initialization
        double xn = searchStart;
        double xStar = xn;
        double fxn = computeObjectiveValue(searchStart);

        double  pn, cn , hn, qn, wn, fwn = 0;


        // Keep finding better approximations.
        while (true) {

            // Regula-Falsi Iteration
            cn = getCn(x0, x1, f0, f1);

            // Convergence test
            final double fcn = computeObjectiveValue(cn);
            if(fcn == 0|| FastMath.abs(fcn)<= ftol){
                return cn;
            }

            // check convergence of cn sequence

            if (FastMath.abs(cn - xStar) < FastMath.max(rtol * FastMath.abs(cn), atol)) {
                return cn;
            }


            if(f0 * fcn < 0){

                anBar = x0;
                bnBar = cn;
                fbnBar = fcn;
            }else{
                anBar = cn;
                fanBar = fcn;
                bnBar = x1;
                fbnBar = f1;

            }

            xStar = cn;

            //HEXRF Iteration
            qn = getQn(fxn, f1, f0);
            hn = getHn(x1, x0, fxn, qn);
            pn = getPn(xn, fxn, hn);
            double numerator = qn * FastMath.abs(fxn);
            double denominator = xn * (pn * fxn * fxn + fxn- fcn);
            wn = xn * FastMath.exp(-numerator/denominator);


            if(wn >= anBar && wn <= bnBar) {
                xn = wn;
                fxn = fwn = computeObjectiveValue(wn);
                incrementEvaluationCount();

                if (fanBar * fwn < 0) {
                    x0 = anBar;
                    f0 = fanBar;
                    x1 = wn;
                    f1 = fwn;

                } else {

                    x0 = wn;
                    f0 = fwn;
                    x1 = bnBar;
                    f1 = fbnBar;
                }
            }else{

                x0 = anBar;
                f0 = fanBar;
                x1 = bnBar;
                f1 = fbnBar;
                if(wn < anBar){
                    xn = anBar;
                    fxn = fanBar;
                }else{
                    xn = bnBar;
                    fxn = fbnBar;
                }
            }

            // convergence test
            if( fxn == 0|| FastMath.abs(fxn)<= ftol){

                return xn;
            }

            if (FastMath.abs(x1 - x0) < FastMath.max(rtol * FastMath.abs(x1),
                    atol)) {
               return x1;
            }

            // if xn and xn+1 are close

        }
    }
}
