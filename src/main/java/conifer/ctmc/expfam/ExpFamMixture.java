package conifer.ctmc.expfam;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;
import briefj.Indexer;
import briefj.collections.Counter;

import conifer.ctmc.CTMCParameters;
import conifer.ctmc.RateMatrixToEmissionModel;
import conifer.ctmc.RateMatrixUtils;
import conifer.ctmc.SimpleRateMatrix;
import conifer.ctmc.expfam.features.IdentityBivariate;
import conifer.ctmc.expfam.features.IdentityUnivariate;
import conifer.io.Indexers;
import conifer.io.PhylogeneticObservationFactory;
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
  
  public static ExpFamMixture dnaGTR(int nCategories)
  {
    CTMCExpFam<CTMCState> globalExponentialFamily = CTMCExpFam.createModelWithFullSupport(simpleDNAStateIndexer(nCategories), true);
    globalExponentialFamily.extractReversibleBivariateFeatures(Collections.singleton(new IdentityBivariate<CTMCState>()));
    globalExponentialFamily.extractUnivariateFeatures(Collections.singleton(new IdentityUnivariate<CTMCState>()));
    
    ExpFamParameters params = new ExpFamParameters(globalExponentialFamily);
    RateMatrixToEmissionModel emissionModel = null;
    Indexer<String> dnaIndexer = Indexers.dnaIndexer();
    CTMCStateSpace stateSpace = new CTMCStateSpace(dnaIndexer, dnaIndexer, nCategories);
    ExpFamMixture mixture = new ExpFamMixture(params, emissionModel, stateSpace);
    return mixture;
  }
  
  // TODO: add back indexers in CTMCParameters, create checks to make
  // sure they correspond to the ones used to load the data
  
//  private static class PartitionCache
//  {
//    private final List<Double> logPriors;
//    private final List<CTMCParameters> rateMatrices;
//  }
  
  public ExpFamMixture(
      ExpFamParameters parameters,
      RateMatrixToEmissionModel emissionModel, 
      CTMCStateSpace stateSpace)
  {
    this.parameters = parameters;
    this.emissionModel = emissionModel;
    this.stateSpace = stateSpace;
  }

  private CTMCExpFam<CTMCState>.LearnedReversibleModel getModel()
  {
    // TODO: cache
    CTMCExpFam<CTMCState>.LearnedReversibleModel result = parameters.globalExponentialFamily.reversibleModelWithParameters(parameters.getVector());
    return result;
  }
  
  @Override
  public CTMCParameters getRateMatrix(int categoryIndex)
  {
    final int nLatentStates = stateSpace.latentIndexer.size();
    double [][] rateMatrix = new double[nLatentStates][nLatentStates];
    
    
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
    
    SimpleRateMatrix result = new SimpleRateMatrix(rateMatrix, emissionModel);
    
    // TODO: cache!
    
    // TODO: will need some special treatment for invariant case
    
    return result;
  }

  @Override
  public List<Double> getLogPriorProbabilities()
  {
    CTMCExpFam<CTMCState>.LearnedReversibleModel model = getModel();
    
    // TODO: rename all instance of statio to quasi-station in exp fam stuff
    
    List<Double> result = Lists.newArrayList();
    
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
      result.add(Math.log(sum));
    }
    
    // TODO: cache!
    
    return result;
  }

}
