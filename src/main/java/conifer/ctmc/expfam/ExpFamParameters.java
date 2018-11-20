package conifer.ctmc.expfam;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import com.google.common.collect.Lists;

import blang.core.WritableRealVar;
import conifer.ctmc.RateMatrixUtils;

import briefj.Indexer;
import briefj.collections.Counter;


public class ExpFamParameters 
{
    public final CTMCExpFam<CTMCState> globalExponentialFamily;
    public Counter<Object> weights;
    private final CTMCStateSpace stateSpace;
    
    public WritableRealVar getRealVar(Object key) {
      final Counter<Object> weights = this.weights;
      return new WritableRealVar() {
        @Override
        public double doubleValue() {
          return weights.getCount(key);
        }

        @Override
        public void set(double value) {
          weights.setCount(key, value);
        }
      };
    }

    private CTMCExpFam<CTMCState>.LearnedReversibleModel _cachedModel = null;
    public CTMCExpFam<CTMCState>.LearnedReversibleModel getModel()
    {
        if (_cachedModel != null) return _cachedModel;
        _cachedModel = globalExponentialFamily.reversibleModelWithParameters(getVector());
        return _cachedModel;
    }

    private void clearCache()
    {
        this._cachedModel = null;
        this._rateCache = null;
        this._cachedLogPriorPrs = null;
    }

    private double[][][] _rateCache = null;
    public double[][] getRateMatrix(int categoryIndex)
    {
        if (_rateCache != null && _rateCache[categoryIndex] != null)
            return _rateCache[categoryIndex];
        if (_rateCache == null)
            _rateCache = new double[stateSpace.nCategories][][];
        final int nLatentStates = stateSpace.latentIndexer.size();
        double [][] rateMatrix = new double[nLatentStates][nLatentStates];
        _rateCache[categoryIndex] = rateMatrix;

        CTMCExpFam<CTMCState>.LearnedReversibleModel model = getModel();
        for (int latentStateIndex = 0; latentStateIndex < nLatentStates; latentStateIndex++)
        {
            Object latentState = stateSpace.latentIndexer.i2o(latentStateIndex);
            CTMCState source = new CTMCState(categoryIndex, latentState, stateSpace.currentPartition);
            Counter<CTMCState> rates = model.getRates(source);

            for (int latentStateIndex2 = 0; latentStateIndex2 < nLatentStates; latentStateIndex2++)
                if (latentStateIndex != latentStateIndex2)
                {
                    Object latentState2 = stateSpace.latentIndexer.i2o(latentStateIndex2);
                    CTMCState dest = new CTMCState(categoryIndex, latentState2, stateSpace.currentPartition);
                    double rate = rates.getCount(dest);
                    rateMatrix[latentStateIndex][latentStateIndex2] = rate;
                }
        }
        RateMatrixUtils.fillRateMatrixDiagonalEntries(rateMatrix);
        return rateMatrix;
    }

    private List<Double> _cachedLogPriorPrs = null;
    public List<Double> getLogPriorProbabilities()
    {
        if (_cachedLogPriorPrs != null)
            return _cachedLogPriorPrs;
        CTMCExpFam<CTMCState>.LearnedReversibleModel model = getModel();

        _cachedLogPriorPrs = Lists.newArrayList();

        Counter<CTMCState> quasiStationary = model.getStationaryDistribution();
        final int nLatentStates = stateSpace.latentIndexer.size();
        for (int category = 0; category < stateSpace.nCategories; category++)
        {
            double sum = 0.0;

            for (int latentStateIndex = 0; latentStateIndex < nLatentStates; latentStateIndex++)
            {
                Object latentState = stateSpace.latentIndexer.i2o(latentStateIndex);
                sum += quasiStationary.getCount(new CTMCState(category, latentState, stateSpace.currentPartition));
            }
            _cachedLogPriorPrs.add(Math.log(sum));
        }

        return _cachedLogPriorPrs;
    }


    public ExpFamParameters(CTMCExpFam<CTMCState> globalExponentialFamily, CTMCStateSpace space)
    {
        this(globalExponentialFamily, new Counter<Object>(), space);
    }

    public ExpFamParameters(CTMCExpFam<CTMCState> globalExponentialFamily,
                            Counter<Object> weights, CTMCStateSpace space)
    {
        this.globalExponentialFamily = globalExponentialFamily;
        this.weights = weights;
        this.stateSpace = space;
    }

    public double [] getVector()
    {
        return globalExponentialFamily.convertFeatureCounter(this.weights);
    }

    public void setVector(double [] values)
    {
        clearCache();
        Indexer<Object> fIndexer = globalExponentialFamily.featuresIndexer;
        for (int i = 0; i < values.length; i++){
            weights.setCount(fIndexer.i2o(i), values[i]);
        }

    }

    /**
     * This method will help set the values of weights based on the hashmap
     * featureWeights. This method is mainly used when we want to simulate
     * fasta file from a fixed rate matrix given the values of the weights
     * @param featureWeights
     */

    public void setVector(Map<String, Double> featureWeights){
        clearCache();
        Indexer<Object> fIndexer = globalExponentialFamily.featuresIndexer;
        for (int i = 0; i < featureWeights.keySet().size(); i++){
            weights.setCount(fIndexer.i2o(i), featureWeights.get(fIndexer.i2o(i)));
        }
    }


 
    public int getDim()
    {
        return globalExponentialFamily.nFeatures();
    }

    @Override
    /**
     * Row major traverse. ACGT
     */
    public String toString() {
        double[][] rateMatrix = this.getRateMatrix(0);
        return "RateMatrix [rateMatrix=" + Arrays.deepToString(rateMatrix) + "]";
    }
}
