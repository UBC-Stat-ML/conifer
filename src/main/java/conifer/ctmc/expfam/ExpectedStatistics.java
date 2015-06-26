package conifer.ctmc.expfam;
import java.util.Arrays;

import bayonet.distributions.Multinomial;
import conifer.ctmc.RateMatrixUtils;
import conifer.ctmc.RateMtxExpectations;
import conifer.models.RateMatrixMixture;


public class ExpectedStatistics<S>
{
    final double [] holdTimes;
    final double [][] nTrans;
    final double [] nInit;

    public final CTMCExpFam<S> model;

    public ExpectedStatistics(CTMCExpFam<S> model)
    {
        this.model = model;
        this.holdTimes = new double[model.nStates];
        this.nInit = new double[model.nStates];
        this.nTrans = new double[model.nStates][];
        for (int state = 0; state < model.nStates; state++)
            nTrans[state] = new double[model.supports[state].length];
    }

//  public void addFromMarginalizedData(EndPointDataset<S> data, double [][] rateMatrix)
//  {
//    Counter<S> initialCounts = data.initialCounts;
//    for (S state : initialCounts.keySet())
//      addInitialValue(state, initialCounts.getCount(state));
//    for (double len : data.branchLengths())
//      addMarginalizedPath(data.getEndPointCounter(len), rateMatrix, len);
//  }

    public void addHoldingTime(S state, double time)
    {
        holdTimes[model.stateIndexer.o2i(state)] += time;
    }

    public void addTransition(S state1, S state2, double count)
    {
        int index1 = model.stateIndexer.o2i(state1);
        int index2 = model.stateIndexer.o2i(state2);
        int[] curSupport = model.supports[index1];
        int supportIdx = Arrays.binarySearch(curSupport, index2);
        if (supportIdx < 0)
            throw new RuntimeException("Transition not in the support: " + state1 + " -> " + state2);
        nTrans[index1][supportIdx] += count;
    }

    public void addInitialValue(S state, double count)
    {
        nInit[model.stateIndexer.o2i(state)] += count;
    }


    public void addMarginalizedPath(double [][] marginalCounts, double [][] rateMtx, double T){
        final int dim = rateMtx.length;
        double [][] auxMtx = new double[2*dim][2*dim];
        double [][] simpleExp = RateMatrixUtils.marginalTransitionMtx(rateMtx, T, RateMatrixUtils.MatrixExponentialAlgorithm.DIAGONALIZATION);
        for(int state1=0;state1<dim;state1++){
            int[] curSupport = model.supports[state1];
            for(int state2Idx=0; state2Idx<curSupport.length+1; state2Idx++){
                boolean isHoldTime = state2Idx== curSupport.length;
                int state2 = isHoldTime? state1:curSupport[state2Idx];
                double sum =0.0;
                double [][] current = RateMtxExpectations._expectations(rateMtx, T, state1, state2, simpleExp, auxMtx);
                for (int x = 0; x < model.stateIndexer.size(); x++) {
                    for (int y = 0; y < model.stateIndexer.size(); y++) {
                        sum += marginalCounts[x][y] * current[x][y];
                    }
                }
                if (isHoldTime)
                    holdTimes[state1] += sum;
                else
                    nTrans[state1][state2Idx] += sum;
            }
        }
    }






//  public void addMarginalizedPath(CounterMap<S, S> endPointCounts, double [][] rateMtx, double T)
//  {
//    _addMarginalizedPath(endPointCounts, null, rateMtx, T);
//  }

//   public void addMarginalizedPath(double [][] marginalCounts, double [][] rateMtx, double T)
//  {
//    _addMarginalizedPath(null, marginalCounts, rateMtx, T);
//  }
//  private void _addMarginalizedPath(CounterMap<S, S> endPointCounts, double [][] marginalCounts, double [][] rateMtx, double T)
//  {
//    if (endPointCounts != null && marginalCounts != null) throw new RuntimeException();
//    boolean useCounter = endPointCounts != null;
//    if (useCounter && endPointCounts.totalCount() == 0.0) return;
//    final int dim = rateMtx.length;
//    double [][] auxMtx = new double[2*dim][2*dim];
//    double [][] simpleExp = RateMtxUtils.marginalTransitionMtx(rateMtx, T);//MatrixFunctions.expm(new DoubleMatrix(rateMtx).mul(T)).toArray2();
//    for (int state1 = 0; state1 < dim; state1++)
//    {
//      int[] curSupport = model.supports[state1];
//      for (int state2Idx = 0; state2Idx < curSupport.length + 1; state2Idx++)
//      {
//        boolean isHoldTime = state2Idx == curSupport.length;
//        int state2 = isHoldTime ? state1 : curSupport[state2Idx];
//        double sum = 0.0;
//        double [][] current = RateMtxExpectations._expectations(rateMtx, T, state1, state2, simpleExp, auxMtx);
//        if (useCounter)
//          for (S s1 : endPointCounts.keySet())
//          {
//            Counter<S> currentCounter = endPointCounts.getCounter(s1);
//            int i1 = model.stateIndexer.o2i(s1);
//            for (S s2 : currentCounter.keySet())
//            {
//              int i2 = model.stateIndexer.o2i(s2);
//              double currentCount = currentCounter.getCount(s2);
//              sum += currentCount * current[i1][i2];
//            }
//          }
//        else
//          for (int x = 0; x < model.stateIndexer.size(); x++)
//            for (int y = 0; y < model.stateIndexer.size(); y++)
//            {
////              if (Double.isNaN(marginalCounts[x][y]) || Double.isNaN(current[x][y]))
////                System.out.println("" + marginalCounts[x][y] + "," + current[x][y]);
//              sum += marginalCounts[x][y] * current[x][y];
//            }
//        if (isHoldTime)
//          holdTimes[state1] += sum;
//        else
//          nTrans[state1][state2Idx] += sum;
//      }
//    }
//  }
//  
//  public void addCategorySpecificMarginalizedPath(double [][] marginalCounts, double [][] rateMtx, double T, CategoryModel catModel, int cat)
//  {
//    final int nObservations = rateMtx.length;
//    double [][] auxMtx = new double[2*nObservations][2*nObservations];
//    double [][] simpleExp = RateMtxUtils.marginalTransitionMtx(rateMtx, T); //MatrixFunctions.expm(new DoubleMatrix(rateMtx).mul(T)).toArray2();
//    for (int observationIndex1 = 0; observationIndex1 < nObservations; observationIndex1++)
//    {
//      AnnotatedCharacter annChar1 = new AnnotatedCharacter(catModel.observationsIndexer.i2o(observationIndex1), cat);
//      int state1 = catModel.indexer.o2i(annChar1);
//      int[] curSupport = model.supports[state1];
//      for (int stateIdx2 = 0; stateIdx2 < curSupport.length + 1; stateIdx2++)
//      {
//        boolean isHoldTime = stateIdx2 == curSupport.length;
//        int observationIndex2 = isHoldTime ? observationIndex1 : catModel.observationsIndexer.o2i(catModel.indexer.i2o(curSupport[stateIdx2]).observedChar);
//        double sum = 0.0;
//        double [][] current = RateMtxExpectations._expectations(rateMtx, T, observationIndex1, observationIndex2, simpleExp, auxMtx);
//        for (int x = 0; x < nObservations; x++)
//          for (int y = 0; y < nObservations; y++)
//            if (marginalCounts[x][y] != 0.0)
//              sum += marginalCounts[x][y] * current[x][y];
//        if (isHoldTime)
//          holdTimes[state1] += sum;
//        else
//          nTrans[state1][stateIdx2] += sum;
//      }
//    }
//  }

    public double nSeries()
    {
        return Multinomial.getNormalization(nInit);
    }

    public double totalTime()
    {
        return Multinomial.getNormalization(holdTimes);
    }

    @Override
    public String toString()
    {
        return "holdTimes = " + Arrays.toString(holdTimes) + "\nnInitStates = " + Arrays.toString(nInit) + "\nnTransitions =\n" + Arrays.deepToString(nTrans);
    }

//  public void addInitialAndFullyObservedPathStatistics(
//      List<Pair<Integer, Double>> datum)
//  {
//    addInitialValue(model.stateIndexer.i2o(datum.get(0).getFirst()),1);
//    for (int d = 0; d < datum.size(); d++)
//    {
//      S char1 = model.stateIndexer.i2o(datum.get(d).getFirst());
//      addHoldingTime(char1, datum.get(d).getSecond());
//      if (d != datum.size() - 1)
//      {
//        S char2 = model.stateIndexer.i2o(datum.get(d+1).getFirst());
//        addTransition(char1, char2, 1);
//      }
//    }
//  }



//  public static double[][] thetaMLE()
//  {
//    Counter<UnorderedPair<Integer, Integer>> N
//  }
}