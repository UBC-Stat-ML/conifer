package conifer.moves;

import hmc.AHMC;
import hmc.DataStruct;
import hmc.HMC;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.jblas.DoubleMatrix;

import com.google.common.collect.Lists;

import conifer.SingleProteinModel;
import conifer.ctmc.PathStatistics;
import conifer.ctmc.expfam.CTMCExpFam;
import conifer.ctmc.expfam.CTMCState;
import conifer.ctmc.expfam.CTMCStateSpace;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.ExpFamParameters;
import conifer.ctmc.expfam.ExpectedStatistics;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.NodeMove;
import blang.mcmc.SampledVariable;
import briefj.BriefIO;
import briefj.Indexer;
import briefj.opt.Option;
import briefj.run.Results;



public class PhyloHMCMove extends NodeMove
{
  @SampledVariable ExpFamParameters parameters;
  
  @ConnectedFactor UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood;
  @ConnectedFactor IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior;
  
  @Option(gloss="provide epsilon")
  public static Double epsilon = null;
  @Option(gloss="provide L")
  public static Integer L = null;
  
  public int nItersPerPathAuxVar = 1000;

  
  private final PrintWriter detailWriter = BriefIO.output(Results.getFileInResultFolder("HMC.experiment.details.txt"));
  
  @Override
  public void execute(Random rand)
  {
    // hack for now to make this sampled less often
    if (rand.nextInt(10) != 0)
      return;
    
    if (prior.marginalDistributionParameters.mean.getValue() != 0.0)
      throw new RuntimeException();
    final double variance = prior.marginalDistributionParameters.variance.getValue();
    
    List<PathStatistics> pathStatistics = likelihood.evolutionaryModel.samplePosteriorPaths(rand, likelihood.observations, likelihood.tree);
    
    ExpectedStatistics<CTMCState> convertedStat = convert(pathStatistics, parameters, likelihood); 
    CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective objective = parameters.globalExponentialFamily.getExpectedCompleteReversibleObjective(1.0/variance, convertedStat);
    
    double [] initialPoint = parameters.getVector();
    double [] newPoint;
    
    if (hyperParametersInitialized())
    {
      // TODO: if needed, could also do several HMCs accept-reject rounds keeping the 
      // expected stats fixed,
      // but may not be needed (already quite a bit of gains by doing large number of steps (L) 
      // within the doIter() method below
      DataStruct hmcResult = null;
      for (int i = 0; i < nItersPerPathAuxVar; i++)
        hmcResult = HMC.doIter(rand, L, epsilon, i == 0 ? new DoubleMatrix(initialPoint) : hmcResult.next_q, objective, objective);
      newPoint = hmcResult.next_q.data;
    }
    else
    {
      AHMC ahmc = AHMC.initializeAHMCWithLBFGS(10000, 1000, objective, objective, initialPoint.length,1000);
      newPoint = ahmc.sample(rand).data;
      epsilon = ahmc.getEpsilon();
      L = ahmc.getL();
      logToFile("optimal epsilon" +""+ epsilon);
      logToFile("optimal L" + ""+ L);
    }
    
    parameters.setVector(newPoint);
  }

  private boolean hyperParametersInitialized()
  {
    return epsilon != null;
  }

  public static ExpectedStatistics<CTMCState> convert(
      List<PathStatistics> pathStatistics,
      ExpFamParameters parameters,
      UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood)
  {
    ExpectedStatistics<CTMCState> result = new ExpectedStatistics<CTMCState>(parameters.globalExponentialFamily);
    CTMCStateSpace space = likelihood.evolutionaryModel.rateMatrixMixture.stateSpace;
    
    for (int category = 0; category < pathStatistics.size(); category++)
    {
      PathStatistics currentStat = pathStatistics.get(category);
      List<CTMCState> states = states(category, space);
      
      for (int s0 = 0; s0 < states.size(); s0++)
      {
        CTMCState state0 = states.get(s0);
        
        // holding times
        result.addHoldingTime(state0, currentStat.getSojournTime(s0));
        
        // initials
        result.addInitialValue(state0, currentStat.getInitialCount(s0));
      
        // transitions
        for (int s1 = 0; s1 < states.size(); s1++)
          if (s0 != s1)
            result.addTransition(state0, states.get(s1), currentStat.getTransitionCount(s0, s1));
      }
      
    }
    
    return result;
  }

  public static List<CTMCState> states(int category, CTMCStateSpace space)
  {
    List<CTMCState> result = Lists.newArrayList();
    
    Object partitionId = space.currentPartition;
    Indexer<?> latentIndexer = space.latentIndexer;
    
    for (int i = 0; i < latentIndexer.size(); i++)
      result.add(new CTMCState(category, latentIndexer.i2o(i), partitionId));
    
    return result;
  }
  
  public void logToFile(String someline) {
    this.detailWriter.println(someline);
    this.detailWriter.flush();
  }

}
