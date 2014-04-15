package conifer.moves;

import hmc.AHMC;
import hmc.DataStruct;
import hmc.HMC;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.jblas.DoubleMatrix;

import com.google.common.collect.Lists;

import conifer.UnrootedTree;
import conifer.ctmc.PathStatistics;
import conifer.ctmc.expfam.CTMCExpFam;
import conifer.ctmc.expfam.CTMCState;
import conifer.ctmc.expfam.CTMCStateSpace;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.ExpFamParameters;
import conifer.ctmc.expfam.ExpectedStatistics;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;

import bayonet.distributions.Gamma;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.NodeMove;
import blang.mcmc.SampledVariable;
import briefj.Indexer;



public class AllBranchesScaling extends NodeMove
{
  @SampledVariable UnrootedTree tree;
  
  @ConnectedFactor UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood;
  @ConnectedFactor NonClockTreePrior<Gamma.Parameters> prior;
  

  @Override
  public void execute(Random rand)
  {
    // hack for now to make this sampled less often
    if (rand.nextInt(10) != 0)
      return;
    
    throw new RuntimeException();
    
//    if (prior.marginalDistributionParameters.mean.getValue() != 0.0)
//      throw new RuntimeException();
//    final double variance = prior.marginalDistributionParameters.variance.getValue();
//    
//    List<PathStatistics> pathStatistics = likelihood.evolutionaryModel.samplePosteriorPaths(rand, likelihood.observations, likelihood.tree);
//    
//    ExpectedStatistics<CTMCState> convertedStat = convert(pathStatistics); 
//    CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective objective = parameters.globalExponentialFamily.getExpectedCompleteReversibleObjective(1.0/variance, convertedStat);
//    
//    double [] initialPoint = parameters.getVector();
//    double [] newPoint;
//    
//    if (hyperParametersInitialized())
//    {
//      // TODO: if needed, could also do several HMCs accept-reject rounds keeping the 
//      // expected stats fixed,
//      // but may not be needed (already quite a bit of gains by doing large number of steps (L) 
//      // within the doIter() method below
//      System.out.println("started HMC");
//      DataStruct hmcResult = null;
//      for (int i = 0; i < 1000; i++)
//      {
//        hmcResult = HMC.doIter(rand, L, epsilon, i == 0 ? new DoubleMatrix(initialPoint) : hmcResult.next_q, objective, objective);
//        System.out.println(hmcResult.energy + "\t" + Arrays.toString(hmcResult.next_q.data));
//      }
//      System.out.println("HMC completed");
//      newPoint = hmcResult.next_q.data;
//    }
//    else
//    {
//      AHMC ahmc = AHMC.initializeAHMCWithLBFGS(3000, 1000, objective, objective, initialPoint.length);
//      newPoint = ahmc.sample(rand).data;
//      epsilon = ahmc.getEpsilon();
//      L = ahmc.getL();
//    }
//    
//    parameters.setVector(newPoint);
  }


}
