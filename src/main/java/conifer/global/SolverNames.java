package conifer.global;

import org.apache.commons.math3.analysis.solvers.BaseAbstractUnivariateSolver;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.analysis.solvers.PegasusSolver;

/**
 * Created by crystal on 2016-10-25.
 */
public enum SolverNames {


    Pegasus{

        public PegasusSolver getSolver(){

            final PegasusSolver solver = new PegasusSolver();
            return solver;
        }




    },

    Brent{

        public BrentSolver getSolver(){

            final BrentSolver solver = new BrentSolver();
            return solver;
        }

    };

    public abstract BaseAbstractUnivariateSolver getSolver();


}
