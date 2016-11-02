package conifer.global;

import org.apache.commons.math3.analysis.solvers.BaseAbstractUnivariateSolver;
import org.apache.commons.math3.analysis.solvers.PegasusSolver;

/**
 * Created by crystal on 2016-11-01.
 */
public enum SelfImplementedSolverNames {

    Quadratic{

        public QuadraticRegulaFalsiSolver getSolver(){

            final QuadraticRegulaFalsiSolver solver = new QuadraticRegulaFalsiSolver();
            return solver;
        }

    },


    ThirdOrder{

        public ThirdOrderRegulaFalsiSolver getSolver(){

            final ThirdOrderRegulaFalsiSolver solver = new ThirdOrderRegulaFalsiSolver();
            return solver;
        }
    },

    ImprovedPegasus{

        public ImprovedPegasusSolver getSolver(){
            final ImprovedPegasusSolver solver = new ImprovedPegasusSolver();
            return solver;

        }

    };

    public abstract SelfImplementedSolver getSolver();



}
