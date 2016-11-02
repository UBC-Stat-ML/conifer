package conifer.global;


import bayonet.math.EJMLUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BaseAbstractUnivariateSolver;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.analysis.solvers.PegasusSolver;
import org.apache.commons.math3.optim.linear.SolutionCallback;
import org.jblas.DoubleMatrix;

import conifer.global.GlobalRFSampler.CollisionSolver;
import bayonet.math.NumericalUtils;
import bayonet.opt.DifferentiableFunction;
import bayonet.opt.LBFGSMinimizer;

import static conifer.rejfreeutil.StaticUtils.*;

public class MultipleConvexCollisionSolver implements CollisionSolver{

    //private final PegasusSolver solver = new PegasusSolver();
    private final BaseAbstractUnivariateSolver solver;

    private final QuadraticRegulaFalsiSolver quadraticSolver;
    private final ThirdOrderRegulaFalsiSolver thirdOrderRegulaFalsiSolver;
    private final ImprovedPegasusSolver improvedPegasusSolver;

    public static boolean useQuadraticSolver = false;
    public static boolean useThirdOrderSolver = false;
    public static boolean useImprovedPegasusSolver = true;

    public MultipleConvexCollisionSolver(){

        this.solver = SolverNames.Pegasus.getSolver();
        this.quadraticSolver = new QuadraticRegulaFalsiSolver();
        this.thirdOrderRegulaFalsiSolver = new ThirdOrderRegulaFalsiSolver();
        this.improvedPegasusSolver = new ImprovedPegasusSolver();
    }

    public MultipleConvexCollisionSolver(SolverNames selectedSolver){

        this.solver = selectedSolver.getSolver();
        this.quadraticSolver = new QuadraticRegulaFalsiSolver();
        this.thirdOrderRegulaFalsiSolver = new ThirdOrderRegulaFalsiSolver();
        this.improvedPegasusSolver = new ImprovedPegasusSolver();
    }

    public double collisionTime(final DoubleMatrix initialPoint, final DoubleMatrix velocity, DifferentiableFunction energy, final double exponential) {
        // go to minimum energy for free
        final DoubleMatrix directionalMin = lineMinimize(initialPoint, velocity, energy);//new DoubleMatrix(searcher.minimize(energy, initialPoint.data, velocity.data));

        // When this is used in the local version of the sampler, there can be situations
        // where one of the factor is improper, causing some trajectory to not incur collision
        // with respect to one of the factors. As long as the product of all the factors in
        // proper, the other factors will ensure that all the variables
        // still get updated infinitely often
        if (directionalMin == null)
            return Double.POSITIVE_INFINITY;

        final double time1 = time(initialPoint, directionalMin, velocity);

        // keep moving until an exponentially distributed amount of energy is exhausted

        final double initialEnergy = energy.valueAt(directionalMin.data);
        final UnivariateFunction lineSolvingFunction = new UnivariateFunction() {
            @Override
            public double value(final double time) {
                final DoubleMatrix candidatePosition = position(directionalMin, velocity, time);
                final double candidateEnergy = energy.valueAt(candidatePosition.data);
                final double delta = candidateEnergy - initialEnergy;
                if (delta < -NumericalUtils.THRESHOLD)
                    System.err.println("Did not expect negative delta for convex objective. " +
                            "Delta=" + delta + ", time=" + time);
                return exponential - delta;
            }
        };
        final double upperBound = findUpperBound(lineSolvingFunction);
        final int maxEval = 100;
        if(useQuadraticSolver){
//            final double time3 = solver.solve(maxEval, lineSolvingFunction, 0.0, upperBound);
            double time2 = quadraticSolver.solve(maxEval,lineSolvingFunction, 0.0, upperBound);

//            if(Math.abs(time2-time3)>1e-6){
//                System.out.println(Double.toString(time3));
//                System.out.println(Double.toString(time2));
////                //throw new RuntimeException("Two solvers get different solution");
//           }
            return time1 + time2;

        }

        else if(useThirdOrderSolver){
            //final double time3 = solver.solve(maxEval, lineSolvingFunction, 0.0, upperBound);
            double time2 = thirdOrderRegulaFalsiSolver.solve(maxEval,lineSolvingFunction, 0.0, upperBound);
//            if(Math.abs(time2-time3)>1e-6){
//                System.out.println(Double.toString(time3));
//                System.out.println(Double.toString(time2));
//                //throw new RuntimeException("Two solvers get different solution");
//            }
            return time1 + time2;

        }
        else if(useImprovedPegasusSolver){

            final double time3 = solver.solve(maxEval, lineSolvingFunction, 0.0, upperBound);
            double time2 = improvedPegasusSolver.solve(3 * maxEval,lineSolvingFunction, 0.0, upperBound);
            if(Math.abs(time2-time3)>1e-6){
                System.out.println(Double.toString(time3));
                System.out.println(Double.toString(time2));
                //throw new RuntimeException("Two solvers get different solution");
            }
            return time1 + time2;



        }

        else{
             double time2 = solver.solve(maxEval, lineSolvingFunction, 0.0, upperBound);
             return time1 + time2;
            }



    }

    private static DoubleMatrix lineMinimize(
            final DoubleMatrix initialPoint,
            final DoubleMatrix velocity,
            final DifferentiableFunction energy) {
        DifferentiableFunction lineRestricted = new DifferentiableFunction() {

            @Override
            public double valueAt(double[] _time) {
                double time = _time[0];
                double[] position = position(initialPoint, velocity, time).data;
                return energy.valueAt(position);
            }

            @Override
            public int dimension() {
                return 1;
            }

            @Override
            public double[] derivativeAt(double[] _time) {
                double time = _time[0];
                double[] position = position(initialPoint, velocity, time).data;
                DoubleMatrix fullDerivative = new DoubleMatrix(energy.derivativeAt(position));
                double directionalDeriv = fullDerivative.dot(velocity);
                return new double[]{directionalDeriv};
            }
        };

        double minTime = new LBFGSMinimizer().minimize(lineRestricted, new double[]{0}, 1e-10)[0];

        if (minTime < 0.0)
            minTime = 0.0;

        double minValue = lineRestricted.valueAt(new double[]{minTime});
        double valuePlusDelta = lineRestricted.valueAt(new double[]{minTime + DELTA});
        if (valuePlusDelta < minValue)
            return null;

        return position(initialPoint, velocity, minTime);
    }

    private static final double DELTA = 1.0;

    private static double time(DoubleMatrix initialPos, DoubleMatrix finalPosition, DoubleMatrix velocity) {
        final double
                xInit = initialPos.get(0),
                xFinal = finalPosition.get(0),
                v = velocity.get(0);
        return (xFinal - xInit) / v;
    }

    private static double findUpperBound(UnivariateFunction lineSolvingFunction) {
        double result = 1.0;
        final int maxNIterations = Double.MAX_EXPONENT - 1;
        for (int i = 0; i < maxNIterations; i++) {
            if (lineSolvingFunction.value(result) < 0.0)
                return result;
            else
                result *= 2.0;
        }
        throw new RuntimeException();
    }
}
