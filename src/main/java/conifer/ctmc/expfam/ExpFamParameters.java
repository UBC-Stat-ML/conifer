package conifer.ctmc.expfam;

import conifer.moves.RealVectorMHProposal;
import blang.annotations.Processors;
import blang.annotations.Samplers;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.variables.RealVectorInterface;
import briefj.Indexer;
import briefj.collections.Counter;


@Samplers(RealVectorMHProposal.class)
@Processors(IIDRealVectorGenerativeFactor.VectorNormProcessor.class)
public class ExpFamParameters implements RealVectorInterface
{
  public final CTMCExpFam<CTMCState> globalExponentialFamily;
  private Counter<Object> weights;
  

  
//  private Map<Pair<Integer,Object>, CTMCParameters> rateMatrices = null;
//  private Map<Object,List<Double>> logPriorPrs = null;
  
//  private final Map cache = Maps.newHashMap();
  
  // TODO: need some kind of caching construct?
  
//  private List<double[][]> matrices;
//  private List<Double> logWeights;
  
  public ExpFamParameters(CTMCExpFam<CTMCState> globalExponentialFamily)
  {
    this(globalExponentialFamily, new Counter<Object>());
  }
  
  public ExpFamParameters(CTMCExpFam<CTMCState> globalExponentialFamily,
      Counter<Object> weights)
  {
    this.globalExponentialFamily = globalExponentialFamily;
    this.weights = weights;
  }

  public double [] getVector()
  {
    return globalExponentialFamily.convertFeatureCounter(this.weights);
  }
  
  public void setVector(double [] values)
  {
    Indexer<Object> fIndexer = globalExponentialFamily.featuresIndexer;
    for (int i = 0; i < values.length; i++)
      weights.setCount(fIndexer.i2o(i), values[i]);
  }
//  
//  private void clearCache()
//  {
//    cache.clear();
//  }
  


//  public CTMCParameters getRateMatrix(int categoryIndex, ExpFamStateSpace stateSpace)
//  {
//    Object currentPartition = stateSpace.getPartition()
//    ensureCache(currentPartition);
//    return cache.get(currentPartition).rateMatrices.get(categoryIndex); 
//  }
//
//  public List<Double> getLogPriorProbabilities(Object currentPartition)
//  {
//    ensureCache(currentPartition);
//    return cache.get(currentPartition).logPriors;
//  }
//
//  private void ensureCache(Object partition)
//  {
//    if (isCacheAvailableFor(partition))
//      return;
//    
//    CTMCExpFam<CTMCState>.LearnedReversibleModel modelWithParams = globalExponentialFamily.reversibleModelWithParameters(weights);
//    
//  }
//  
//  private boolean isCacheAvailableFor(Object partition)
//  {
//    return cache.containsKey(partition);
//  }
  

  
  // weight -> fullMtx(implicitly) -> splitMtx: needs features, support, isNorm
  
  // 
  
  
}
