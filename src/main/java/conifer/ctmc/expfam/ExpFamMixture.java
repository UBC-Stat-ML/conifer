package conifer.ctmc.expfam;

import java.util.Collections;
import java.util.List;

import blang.inits.DesignatedConstructor;
import briefj.Indexer;
import conifer.ctmc.CTMCParameters;
import conifer.ctmc.RateMatrixToEmissionModel;
import conifer.ctmc.SimpleRateMatrix;
import conifer.ctmc.expfam.features.IdentityBivariate;
import conifer.ctmc.expfam.features.IdentityUnivariate;
import conifer.io.Indexers;
import conifer.models.RateMatrixMixture;
import conifer.ctmc.expfam.RateMtxNames;



public class ExpFamMixture implements RateMatrixMixture
{
  
  
  public final ExpFamParameters parameters;
  

  public final RateMatrixToEmissionModel emissionModel;
  
  public final CTMCStateSpace stateSpace;
  
  public static int nFeatures;
  public int nFeatures() { return nFeatures; }
  
  
  
  public static Indexer<CTMCState> simpleDNAStateIndexer(int nCategories, final RateMtxNames selectedRateMtx)
  {
    Indexer<String> indexer= Indexers.modelIndexer(selectedRateMtx);
    Indexer<CTMCState> result = new Indexer<CTMCState>();
    
    for (int cat = 0; cat < nCategories; cat++)
      for (int i = 0; i < indexer.size(); i++)
        result.addToIndex(new CTMCState(cat, indexer.i2o(i), null));
    return result;
  }
  
  public static ExpFamMixture fromSerialized(SerializedExpFamMixture serialized, Indexer<String> observationIndexer, boolean isNormalized)
  {
    CTMCExpFam<CTMCState> globalExponentialFamily = new CTMCExpFam<CTMCState>(serialized.getSupport(), serialized.getCTMCStateIndexer(), isNormalized);
    globalExponentialFamily.extractReversibleBivariateFeatures(Collections.singletonList(serialized.getBivariateFeatureExtractor()));
    globalExponentialFamily.extractUnivariateFeatures(Collections.singletonList(serialized.getUnivariateFeatureExtractor()));
    nFeatures = globalExponentialFamily.nFeatures();
    RateMatrixToEmissionModel emissionModel = serialized.getEmissionModel();
    Indexer<String> latentIndexer = serialized.getLatentIndexer();
    if (emissionModel == null)
      if (!latentIndexer.equals(observationIndexer))
        throw new RuntimeException();
    CTMCStateSpace stateSpace = new CTMCStateSpace(observationIndexer, latentIndexer, serialized.nCategories, SerializedExpFamMixture.partition);
    ExpFamParameters params = new ExpFamParameters(globalExponentialFamily, stateSpace);
    ExpFamMixture mixture = new ExpFamMixture(params, emissionModel, stateSpace);
    return mixture;
  }
  
  
  
  public static ExpFamMixture rateMtxModel(RateMtxNames selectedRateMtx, boolean isNormalized)
  {
   return fromSerialized(SerializedExpFamMixture.rateMtxModel(selectedRateMtx), Indexers.modelIndexer(selectedRateMtx), isNormalized);
     }
  
 
  public static ExpFamMixture randomGTR(RateMtxNames selectedRateMtx)
  {
    final int nCat = 1;
    CTMCExpFam<CTMCState> globalExponentialFamily = CTMCExpFam.createModelWithFullSupport(simpleDNAStateIndexer(nCat,selectedRateMtx), true);
    globalExponentialFamily.extractReversibleBivariateFeatures(Collections.singleton(new IdentityBivariate<CTMCState>()));
    globalExponentialFamily.extractUnivariateFeatures(Collections.singleton(new IdentityUnivariate<CTMCState>()));
    
    RateMatrixToEmissionModel emissionModel = null;
    Indexer<String> dnaIndexer = Indexers.modelIndexer(selectedRateMtx);
    CTMCStateSpace stateSpace = new CTMCStateSpace(dnaIndexer, dnaIndexer, nCat);
    ExpFamParameters params = new ExpFamParameters(globalExponentialFamily, stateSpace);
    ExpFamMixture mixture = new ExpFamMixture(params, emissionModel, stateSpace);
    return mixture;
  }
  

  public ExpFamMixture(
      ExpFamParameters parameters,
      RateMatrixToEmissionModel emissionModel, 
      CTMCStateSpace stateSpace)
  {
    this.parameters = parameters;
    this.emissionModel = emissionModel;
    this.stateSpace = stateSpace;
  }

  @Override
  public CTMCParameters getRateMatrix(int categoryIndex)
  {
    SimpleRateMatrix result = new SimpleRateMatrix(parameters.getRateMatrix(categoryIndex), emissionModel);
    // TODO: will need some special treatment for invariant case
    return result;
  }

  @Override
  public List<Double> getLogPriorProbabilities()
  {
    return parameters.getLogPriorProbabilities();
  }

}
