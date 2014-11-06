package conifer.ctmc.expfam;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.apache.commons.lang3.tuple.Pair;
import org.jblas.DoubleMatrix;
import org.jgrapht.Graphs;
import org.jgrapht.UndirectedGraph;

import utils.MultiVariateObj;
import utils.Objective;

import bayonet.distributions.Multinomial;
import bayonet.graphs.GraphUtils;
import bayonet.math.SparseVector;
import bayonet.opt.DifferentiableFunction;
import bayonet.opt.LBFGSMinimizer;
import bayonet.opt.OptimizationOptions;
import briefj.Indexer;
import briefj.collections.Counter;
import conifer.ctmc.RateMatrixUtils;



// TODO: make the base measure be a Counter; and center on e.g. multicategory gamma matrix
public class CTMCExpFam<S>
{
  public final Indexer<S> stateIndexer;
  final int [][] supports; // S index -> list of S indices
  private final SparseVector [][] bivariateFeatures; // S index -> index in support
  private final SparseVector [] univariateFeatures;
  final int nStates;
  private int nFeatures;
  public final Indexer<Object> featuresIndexer = new Indexer<Object>();
  public final  boolean  isNormalized;
  public int nFeatures() { return nFeatures; }

  public CTMCExpFam(
      UndirectedGraph<S,?> support, 
      Indexer<S> indexer,
      boolean isNormalized)
  {
    // create support info
    this.stateIndexer = indexer; //= new Indexer<S>(support.vertexSet());
    this.nStates = this.stateIndexer.size();
    this.supports = new int[nStates][];
    this.bivariateFeatures = new SparseVector[nStates][];
    this.univariateFeatures = new SparseVector[nStates];
    this.isNormalized = isNormalized;
    for (int state = 0; state < nStates; state++)
    {
      S current = stateIndexer.i2o(state);
      Collection<S> nbhr = Graphs.neighborListOf(support, current);
      this.supports[state] = new int[nbhr.size()];
      int i = 0;
      for (S item : nbhr)
      {
        if (item == current)
          throw new RuntimeException();
        this.supports[state][i++] = stateIndexer.o2i(item);
      }
      Arrays.sort(this.supports[state]);
      this.bivariateFeatures[state] = new SparseVector[nbhr.size()];
    }
  }

  public static <S> CTMCExpFam<S> createModelWithFullSupport(Indexer<S> indexer, boolean isNormalized)
  {
    return new CTMCExpFam<S>(GraphUtils.completeGraph(indexer.objects()), indexer, isNormalized);
  }

  public void extractUnivariateFeatures(Collection<? extends UnivariateFeatureExtractor<S>> univariateFeatureExtractors)
  {
    _extractFeatures(true,  featuresIndexer, null, univariateFeatureExtractors);
  }
  /**
   * Note that the order of the two arguments in transitionFeatureExtractor calls is made canonical automatically
   * (hence the reversible part of the name of this method)
   */
  public void extractReversibleBivariateFeatures(Collection<? extends BivariateFeatureExtractor<S>> bivariateFeatureExtractors)
  {
    _extractFeatures(false, featuresIndexer, bivariateFeatureExtractors, null);
  }

  public LearnedReversibleModel fitReversibleModel(
      OptimizationOptions optimizationOptions,
      ExpectedStatistics<S> currentStats,
      double [] warmStart)
  {
    checkFeaturesInitialized();
    ExpectedCompleteReversibleObjective expectedCompleteObjective =  getExpectedCompleteReversibleObjective(optimizationOptions.regularizationStrength, currentStats);
    LBFGSMinimizer minimizer = new LBFGSMinimizer(optimizationOptions.maxIterations);
    if (warmStart == null)
      warmStart = new double[nFeatures];
    double [] w = minimizer.minimize(expectedCompleteObjective, warmStart, optimizationOptions.tolerance);
    return reversibleModelWithParameters(w);
  }

  public LearnedReversibleModel reversibleModelWithParameters(double [] w)
  {
    checkFeaturesInitialized();
    return new LearnedReversibleModel(w, this.isNormalized);
  }

  public LearnedReversibleModel reversibleModelWithParameters(Counter<Object> counter)
  {
    return reversibleModelWithParameters(convertFeatureCounter(counter));
  }
  
  public double [] convertFeatureCounter(Counter<Object> counter)
  {
    double [] w = new double[nFeatures];
    for (Object o : counter.keySet()) 
      w[featuresIndexer.o2i(o)] = counter.getCount(o);
    return w;
  }

  /**
   * 
   * @param kappa Precision parameter (reciprocal of variance).
   * @param stats
   * @return
   */
  public ExpectedCompleteReversibleObjective getExpectedCompleteReversibleObjective(double kappa, ExpectedStatistics<S> stats)
  {
    checkFeaturesInitialized();
    return new ExpectedCompleteReversibleObjective(kappa, stats);
  }

  public class ExpectedCompleteReversibleObjective 
    implements 
      DifferentiableFunction,
      MultiVariateObj,
      Objective
  {
    private final double kappa; // regularization strength

    private final double [] holdTimes;

    private final double [] nInit;
    private final double nInitStar;

    private final double [][] nTrans; // S index -> index in support
    private final double [] nTransStar;
    private final double nTransStarStar;

    private final double [] fixedDerivative;

    private double lastValue;
    private double[] lastDerivative;
    double[] lastX = null;

    private ExpectedCompleteReversibleObjective(double kappa, ExpectedStatistics<S> stats)
    {
      this.kappa = kappa;
      this.holdTimes = stats.holdTimes;
      this.nInit = stats.nInit;
      this.nTrans = stats.nTrans;

      this.nInitStar = stats.nSeries();
      this.nTransStar = new double[nStates];
      for (int startState = 0; startState < nStates; startState++)
        for (int endStateIdx = 0; endStateIdx < supports[startState].length; endStateIdx++)
        {
          final int endState = supports[startState][endStateIdx];
          nTransStar[endState] += nTrans[startState][endStateIdx];
        }
      this.nTransStarStar = Multinomial.getNormalization(nTransStar); 

      fixedDerivative = _fixedDerivative();
    }

    @Override
    public int dimension()
    {
      return nFeatures;
    }

    public double valueAt(double[] x) {
      ensureCache(x);
      return lastValue;
    }
    public double[] derivativeAt(double[] x) {
      ensureCache(x);
      return lastDerivative;
    }
    private void ensureCache(double[] x) {
      if (requiresUpdate(lastX, x)) {
        Pair<Double, double[]> currentValueAndDerivative = calculate(x);

        lastValue = currentValueAndDerivative.getLeft();
        lastDerivative = currentValueAndDerivative.getRight();
        if (lastX == null)
          lastX = new double[x.length];
        for (int i = 0; i < x.length; i++)
          lastX[i] = x[i];
      }
    }

    private double [] _fixedDerivative()
    {
      double [] result = new double[nFeatures];
      for (int startState = 0; startState < nStates; startState++)
      {
        univariateFeatures[startState].linearIncrement(nInit[startState] + nTransStar[startState], result); // (4) & (7)
        for (int endStateIdx = 0; endStateIdx < supports[startState].length; endStateIdx++)
          bivariateFeatures[startState][endStateIdx].linearIncrement(nTrans[startState][endStateIdx], result); // (6)
      }
      return result;
    }

    private Pair<Double, double[]> calculate(double[] x)
    {
      final double [] gradient = fixedDerivative.clone();
      double value = 0.0;
      LearnedReversibleModel w = new LearnedReversibleModel(x, isNormalized);
      final double [] mStar = new double[nStates];
      double mStarStar = 0.0;

      for (int startState = 0; startState < nStates; startState++)
      {
        final SparseVector curPiFeatures = univariateFeatures[startState];
        curPiFeatures.linearIncrement(-w.pi[startState] * ( nInitStar + nTransStarStar ), gradient); // (5) & (8)
        value += Math.log(w.pi[startState]) * nInit[startState]; // (1)
        final double currentHold = holdTimes[startState];
        final SparseVector qs = w.qs(startState);
        double sumQs = 0.0;
        final int [] curSupports = supports[startState];
        for (int endStateIdx = 0; endStateIdx < curSupports.length; endStateIdx++)
        {
          final double currentQ = qs.values[endStateIdx];
          value += nTrans[startState][endStateIdx] * Math.log(currentQ); // (2)
          sumQs += currentQ;  // (3)
          final double currentM = currentHold * currentQ;
          bivariateFeatures[startState][endStateIdx].linearIncrement(-currentM, gradient); // (9)
          mStar[curSupports[endStateIdx]] += currentM;
          mStarStar += currentM;
        }
        value = value - sumQs * currentHold; // (3) continued
      }

      for (int startState = 0; startState < nStates; startState++)
        univariateFeatures[startState].linearIncrement(w.pi[startState] * mStarStar - mStar[startState], gradient); // (10) & (11)

      if(isNormalized)  // the matrix needs to be normalized
      {
        for (int startState = 0; startState < nStates; startState++)
        {
          final SparseVector qs = w.qs(startState);
          double sumQs = 0.0;
          final int [] curSupports = supports[startState];

          for (int endStateIdx = 0; endStateIdx < curSupports.length; endStateIdx++)
          {
            final double currentQ = qs.values[endStateIdx];
            univariateFeatures[qs.indices[endStateIdx]].linearIncrement(qs.values[endStateIdx]*w.pi[startState]*(nTransStarStar-mStarStar)*(-1),  gradient);
            bivariateFeatures[startState][endStateIdx].linearIncrement(w.pi[startState]*qs.values[endStateIdx]*(nTransStarStar-mStarStar)*(-1),  gradient);
            sumQs+=currentQ;
          }    
          univariateFeatures[startState].linearIncrement(w.pi[startState]*sumQs*(nTransStarStar-mStarStar)*(-1), gradient);
          univariateFeatures[startState].linearIncrement(2*w.pi[startState]*(mStarStar-nTransStarStar)*(-1), gradient);// derivative beta term
        }
      }
      for (int f = 0; f < nFeatures; f++)
      {
        final double curX = x[f];
        gradient[f] = -(gradient[f] - kappa * curX);   // (12)
        value = value - kappa * curX * curX / 2.0;  // (13)
      }
      value = -value;
      return Pair.of(value, gradient);
    }


    private boolean requiresUpdate(double[] lastX, double[] x) {
      if (lastX == null) return true;
      for (int i = 0; i < x.length; i++) {
        if (lastX[i] != x[i])
          return true;
      }
      return false;
    }

    /**
     * Used by AHMC
     */
    @Override
    public double functionValue(DoubleMatrix vec)
    {
      return valueAt(vec.data);
    }

    /**
     * Used by AHMC
     */
    @Override
    public DoubleMatrix mFunctionValue(DoubleMatrix vec)
    {
      return new DoubleMatrix(dimension(),1, derivativeAt(vec.data));
    }
  }

  private void track(Object c) {}
  private void end_track() {}
  private void logs(Object c) {}

  @SuppressWarnings({ "unchecked", "rawtypes" })
  private void _extractFeatures(
      boolean statio, 
      Indexer featureIndexer,
      Collection<? extends BivariateFeatureExtractor<S>> transitionFeatureExtractors,
      Collection<? extends UnivariateFeatureExtractor<S>> stationaryFeatureExtractors)
  {
    track("Extracting features for the " + (statio ? "stationary distribution" : "transitions"));
    Counter counter = new Counter();
    for (int state = 0 ; state < nStates; state++)
      for (int state2Idx = 0; state2Idx < (statio ? 1 : supports[state].length); state2Idx++)
      {
        List<Integer> indices = new ArrayList<Integer>();
        List<Double> values = new ArrayList<Double>();
        S stateObj = stateIndexer.i2o(state);
        S stateObj2 = statio ? null : stateIndexer.i2o(supports[state][state2Idx]);
        track((statio ? "state" : "trans") + 
            "(" + state + (statio ? "" : "," + supports[state][state2Idx]) + ") = " 
            + stateObj + (statio ? "" : " -> " + stateObj2));
        for (Object extractor : (statio ? stationaryFeatureExtractors : transitionFeatureExtractors))
        {
          track("extractor = " + extractor);
          if (statio)
            ((UnivariateFeatureExtractor) extractor).extract(counter, stateObj);
          else
          {
            S _stateObj = state < supports[state][state2Idx] ? stateObj : stateObj2;
            S _stateObj2= state < supports[state][state2Idx] ? stateObj2: stateObj;
            ((BivariateFeatureExtractor) extractor).extract(counter, _stateObj, _stateObj2);
          }
          for (Object feature : counter.keySet())
          {
            if (!featureIndexer.containsObject(feature))
              featureIndexer.addToIndex(feature);
            int fIndex = featureIndexer.o2i(feature);
            double value = counter.getCount(feature);
            logs("feature(" + fIndex + ") = " + feature + " [value=" + value + "]");
            indices.add(fIndex);
            values.add(value);
          }
          end_track();
          counter.clear();
        }
        logs("nFeatures = " + indices.size());
        if (statio)
          univariateFeatures[state] = new SparseVector(indices, values);
        else
          bivariateFeatures[state][state2Idx] = new SparseVector(indices, values);
        end_track();
      }
    end_track();
    nFeatures = featureIndexer.size();
  }

  public class LearnedReversibleModel
  {
    public final double [] weights;
    public final double [] pi;
    public final double normalization;
    
    private LearnedReversibleModel(double [] w, boolean isNormalized)
    {
      this.weights = w;
      this.pi = _buildPi();
      this.normalization = isNormalized ? normalization(w) : 1.0;
    }

    private double [] _buildPi()
    {
      double [] pi = new double[nStates];
      for (int i = 0; i < nStates; i++)
        pi[i] = univariateFeatures[i].dotProduct(weights);
      Multinomial.expNormalize(pi);
      return pi;
    }

    private double normalization(double [] x)
    {
      if (!isNormalized)
        throw new RuntimeException();
      LearnedReversibleModel w = new LearnedReversibleModel(x, false);
      double betainv = 0;
      for (int startState=0; startState < nStates; startState++)
      {
        final int [] curSupports = supports[startState];
        final SparseVector qs = w.qs(startState);
        double sumQs = 0.0;
        for (int endStateIdx=0; endStateIdx < curSupports.length; endStateIdx++)
          sumQs += qs.values[endStateIdx];
        betainv = betainv + w.pi[startState]*sumQs;
      }
      return 1.0/betainv;
    }

    private SparseVector qs(int startState)
    {
      int [] support = supports[startState];
      final int supportSize = support.length;
      double [] values = new double[supportSize];

      for (int i = 0; i < supportSize; i++)
        values[i] = Math.exp(bivariateFeatures[startState][i].dotProduct(weights)) * pi[support[i]]*normalization;
      return new SparseVector(support, values);  
    }

    public Counter<S> getRates(S source)
    {
      int s = stateIndexer.o2i(source);
      SparseVector qs = qs(s);
      int [] support = supports[s];
      Counter<S> result = new Counter<S>();
      for (int j = 0; j < support.length; j++)
        result.setCount(stateIndexer.i2o(support[j]), qs.values[j]);
      return result;
    }    
    
    public Counter<S> getStationaryDistribution()
    {
      Counter<S> result = new Counter<S>();
      for (int i = 0; i < pi.length; i++)
        result.setCount(stateIndexer.i2o(i), pi[i]);
      return result;
    }

    public double[][] getRateMatrix()
    {
      double [][] result = new double[nStates][nStates];
      for (int s1 = 0 ; s1 < nStates; s1++)
      {
        SparseVector qs = qs(s1);
        int [] support = supports[s1];
        for (int j = 0; j < support.length; j++)
          result[s1][support[j]] = qs.values[j];
      }

      RateMatrixUtils.fillRateMatrixDiagonalEntries(result);
      return result;
    }
    public Counter<String> getWeights()
    {
      Counter<String> result = new Counter<String>();

      for (int i = 0; i < nFeatures; i++)
        result.setCount(featuresIndexer.i2o(i).toString(), weights[i]);

      return result;
    }
  }

  private void checkFeaturesInitialized()
  {
    if (featuresIndexer.size()== 0)
      throw new RuntimeException();
  }
}
