package conifer.ctmc;

/**
 * Using and augmented matrix with a Pade+scaling approach (faster for big matrices and does not require diagonalizable matrices)
 * You can consult the following paper  "Comparison of methods for calculating conditional expectations of sufficient statistics for continuous time Markov chains"
 * by Tataru (2011) to understand how this implementation works in order to get those sufficient statistics
 * We are implementing algorithm 3 "EXPM" of this paper
 * This is a revised class based on RateMtxExpectations.java from the old conifer in the svn of Alex's code
 * @author zhaott0416@gmail.com
 */
public class RateMtxExpectations
{
    /**
     * for (a != b):
     * A[i][j][a][b] = E[N(a->b)|X_0=i, X_T=j]
     * for (a == b):
     * A[i][j][a][a] = E[T(a)|X_0=i, X_T=j]
     * where
     * N(a->b) is the number of transitions from a to b in the interval [0,T] and
     * T(a) is the time spent at time in state a in the interval [0,T]
     * @param T
     * @param rateMtx
     * @return sufficient statistics of the the holding time and transitions
     */
    public static double [][][][] expectations(double [][] rateMtx, double T)
    {
        final int n = rateMtx.length;
        double [][] simpleExp = RateMatrixUtils.marginalTransitionMtx(rateMtx, T);
        double [][][][] result = new double[n][n][n][n];
        double [][] emptyMtx = new double[2*n][2*n];
        for (int state1 = 0; state1 < n; state1++)
            for (int state2 = 0; state2 < n; state2++)
            {
                double [][] current = _expectations(rateMtx, T, state1, state2, simpleExp, emptyMtx);
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        result[i][j][state1][state2] = current[i][j];
            }
        return result;
    }

    public static double[][] expectations(double [][] marginalCounts, double [][] rateMtx, double T)
    {
        final int dim = rateMtx.length;
        double [][] result = new double[dim][dim];
        double [][] auxMtx = new double[2*dim][2*dim];
        double [][] simpleExp = RateMatrixUtils.marginalTransitionMtx(rateMtx, T);//MatrixFunctions.expm(new DoubleMatrix(rateMtx).mul(T)).toArray2();
        for (int state1 = 0; state1 < dim; state1++)
            for (int state2 = 0; state2 < dim; state2++)
            {
                double [][] current = _expectations(rateMtx, T, state1, state2, simpleExp, auxMtx);
                double sum = 0.0;
                for (int i = 0; i < dim; i++)
                    for (int j = 0; j < dim; j++)
                        sum += current[i][j] * marginalCounts[i][j];
                result[state1][state2] = sum;
            }
        return result;
    }

    public static double [][] _expectedWaitingTimes(double [][] rateMtx, double T, int state, double [][] matrixExponential, double[][] emptyMtx)
    {
        return _expectations(rateMtx, T, state, state, matrixExponential, emptyMtx);
    }
    /**
     * If state1 diff than state2, this computes expected counts
     * If state1 equals state2, this computes expected waiting times
     *
     * @param rateMtx
     * @param T
     * @param state1
     * @param state2
     * @param matrixExponential
     * @param emptyMtx
     * @return entry (i,j) of the returned matrix gives the
     */
    public static double [][] _expectations(double [][] rateMtx, double T, int state1, int state2, double [][] matrixExponential, double[][] emptyMtx)
    {
        final int n = rateMtx.length;
        if (rateMtx[0].length != n || matrixExponential.length != n)
            throw new RuntimeException();
        double [][] aux = emptyMtx;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                aux[i][j] = aux[i+n][j+n] = rateMtx[i][j] * T;
        aux[state1][state2+n] = 1.0 * T;

        // if we use RateMatrixUtils.marginalTransitionMtx(double [][] rates, double T) we are using MatrixFunctions.expm().toArray2() in jblas
        double [][] exponentiatedAux = RateMatrixUtils.marginalTransitionMtx(aux, 1.0);

        double [][] result = new double[n][n];

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                result[i][j] = exponentiatedAux[i][n+j] / matrixExponential[i][j] * (state1 == state2 ? 1.0 : rateMtx[state1][state2]);

        aux[state1][state2+n] = 0.0;

        return result;
    }

}
