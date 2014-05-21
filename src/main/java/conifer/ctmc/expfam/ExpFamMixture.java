package conifer.ctmc.expfam;

import java.util.Collections;
import java.util.List;

import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;
import briefj.Indexer;
import conifer.ctmc.CTMCParameters;
import conifer.ctmc.RateMatrixToEmissionModel;
import conifer.ctmc.SimpleRateMatrix;
import conifer.ctmc.expfam.features.IdentityBivariate;
import conifer.ctmc.expfam.features.IdentityUnivariate;
import conifer.io.Indexers;
import conifer.models.RateMatrixMixture;



public class ExpFamMixture implements RateMatrixMixture
{
  
  @FactorArgument
  public final ExpFamParameters parameters;
  
  @FactorComponent
  public final RateMatrixToEmissionModel emissionModel;
  
  public final CTMCStateSpace stateSpace;
  
  public static Indexer<CTMCState> simpleDNAStateIndexer(int nCategories)
  {
    Indexer<String> indexer = Indexers.dnaIndexer();
    Indexer<CTMCState> result = new Indexer<CTMCState>();
    for (int cat = 0; cat < nCategories; cat++)
      for (int i = 0; i < indexer.size(); i++)
        result.addToIndex(new CTMCState(cat, indexer.i2o(i), null));
    return result;
  }
  
  public static Indexer<CTMCState> simpleProteinStateIndexer(int nCategories)
  {
    Indexer<String> indexer = Indexers.proteinIndexer();
    Indexer<CTMCState> result = new Indexer<CTMCState>();
    for (int cat = 0; cat < nCategories; cat++)
      for (int i = 0; i < indexer.size(); i++)
        result.addToIndex(new CTMCState(cat, indexer.i2o(i), null));
    return result;
  }
  
  public static ExpFamMixture fromSerialized(SerializedExpFamMixture serialized, Indexer<String> observationIndexer)
  {
    CTMCExpFam<CTMCState> globalExponentialFamily = new CTMCExpFam<CTMCState>(serialized.getSupport(), serialized.getCTMCStateIndexer(), true);
    globalExponentialFamily.extractReversibleBivariateFeatures(Collections.singletonList(serialized.getBivariateFeatureExtractor()));
    globalExponentialFamily.extractUnivariateFeatures(Collections.singletonList(serialized.getUnivariateFeatureExtractor()));
    
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
  
  public static ExpFamMixture kimura1980()
  {
    return fromSerialized(SerializedExpFamMixture.kimura1980(), Indexers.dnaIndexer());
  }
  
  public static ExpFamMixture Accordance45()
  {
    return fromSerialized(SerializedExpFamMixture.accordance(), Indexers.proteinIndexer());
  }
  
  public static ExpFamMixture dnaGTR()
  {
    final int nCat = 1;
    CTMCExpFam<CTMCState> globalExponentialFamily = CTMCExpFam.createModelWithFullSupport(simpleDNAStateIndexer(nCat), true);
    globalExponentialFamily.extractReversibleBivariateFeatures(Collections.singleton(new IdentityBivariate<CTMCState>()));
    globalExponentialFamily.extractUnivariateFeatures(Collections.singleton(new IdentityUnivariate<CTMCState>()));
    
    RateMatrixToEmissionModel emissionModel = null;
    Indexer<String> dnaIndexer = Indexers.dnaIndexer();
    CTMCStateSpace stateSpace = new CTMCStateSpace(dnaIndexer, dnaIndexer, nCat);
    ExpFamParameters params = new ExpFamParameters(globalExponentialFamily, stateSpace);
    ExpFamMixture mixture = new ExpFamMixture(params, emissionModel, stateSpace);
    return mixture;
  }
  
  public static ExpFamMixture proteinGTR()
  {
    final int nCat = 1;
    CTMCExpFam<CTMCState> globalExponentialFamily = CTMCExpFam.createModelWithFullSupport(simpleProteinStateIndexer(nCat), true);
    globalExponentialFamily.extractReversibleBivariateFeatures(Collections.singleton(new IdentityBivariate<CTMCState>()));
    globalExponentialFamily.extractUnivariateFeatures(Collections.singleton(new IdentityUnivariate<CTMCState>()));
    RateMatrixToEmissionModel emissionModel = null;
    Indexer<String> proteinIndexer = Indexers.proteinIndexer();
    CTMCStateSpace stateSpace = new CTMCStateSpace(proteinIndexer,proteinIndexer, nCat);
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
