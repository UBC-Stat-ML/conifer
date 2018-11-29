package conifer.moves;

import conifer.EvoGLM;
import conifer.TopologyUtils;

import conifer.TreeNode;
import conifer.UnrootedTree;
import hmc.AHMC;
import hmc.DataStruct;
import hmc.HMC;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.time.StopWatch;
import org.jblas.DoubleMatrix;

import com.google.common.collect.Lists;

import bayonet.distributions.Random;
import conifer.ctmc.PathStatistics;
import conifer.ctmc.expfam.CTMCExpFam;
import conifer.ctmc.expfam.CTMCState;
import conifer.ctmc.expfam.CTMCStateSpace;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.ctmc.expfam.ExpFamParameters;
import conifer.ctmc.expfam.ExpectedStatistics;
import conifer.factors.UnrootedTreeLikelihoodUtils;
import conifer.io.TreeObservations;
import conifer.models.MultiCategorySubstitutionModel;
import blang.distributions.Normal;
import blang.core.Factor;
import blang.core.LogScaleFactor;
import blang.core.RealVar;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.SampledVariable;
import blang.mcmc.Sampler;
import blang.mcmc.internals.Callback;
import blang.mcmc.internals.ExponentiatedFactor;
import blang.mcmc.internals.SamplerBuilderContext;
import blang.runtime.internals.objectgraph.Node;
import blang.runtime.internals.objectgraph.StaticUtils;
import briefj.BriefIO;
import briefj.Indexer;
import briefj.run.Results;

public class PhyloHMCMove implements Sampler
{
  @SampledVariable(skipFactorsFromSampledModel = true)
  public EvoGLM fullModel;
  
  @ConnectedFactor // TODO: check they only affect the tree? 
  protected List<LogScaleFactor> numericFactors;
  
  private RealVar annealingParameter;

  public static boolean useAuxiliaryVariable = true;
  
  public static Double epsilon = null;
  
  public static Integer L = null;

  public static Integer sizeAdapt = 500;

  public static int nItersPerPathAuxVar = 1000;
  
  private final PrintWriter detailWriter = BriefIO.output(Results.getFileInResultFolder("HMC.experiment.details.txt"));

  @Override
  public boolean setup(SamplerBuilderContext context) {
    this.annealingParameter = context.getAnnealingParameter();
    return true;
  }

  @Override
  public void execute(Random rand) {
    double variance = fullModel.getVariance().doubleValue();
    
    ExpFamMixture mix = fullModel.getEvolutionaryModel().rateMatrixMixture;
    
	  CTMCExpFam<CTMCState>.ExpectedCompleteReversibleObjective auxObjective = null;
	  CTMCExpFam<CTMCState>.ExpectedReversibleObjectiveUpdateExpectedStat noAuxObjective =null; 
      if(useAuxiliaryVariable){
   	 	List<PathStatistics> pathStatistics = fullModel.getEvolutionaryModel().samplePosteriorPaths(rand, fullModel.getObservations(), fullModel.getTree());
   	 	ExpectedStatistics<CTMCState> convertedStat = convert(pathStatistics, mix.parameters, mix.stateSpace);
        auxObjective = mix.parameters.globalExponentialFamily.getExpectedCompleteReversibleObjective(1.0/variance, convertedStat);
        }else{
       	 ExpectedStatistics<CTMCState> expectedStatistics = fullModel.getEvolutionaryModel().getTotalExpectedStatistics(fullModel.getObservations(), fullModel.getTree(), mix.parameters.globalExponentialFamily);
       	 noAuxObjective = mix.parameters.globalExponentialFamily.getExpectedReversibleObjectiveUpdateExpectedStat(1.0/variance, fullModel.getObservations(),
       	    fullModel.getTree(), fullModel.getEvolutionaryModel(), expectedStatistics, mix.parameters);
       	 }
      
      double [] initialPoint = mix.parameters.getVector();
      double [] newPoint= initialPoint;
      
      if (hyperParametersInitialized())
      {
          // TODO: if needed, could also do several HMCs accept-reject rounds keeping the
          // expected stats fixed,
          // but may not be needed (already quite a bit of gains by doing large number of steps (L)
          // within the doIter() method below
          DataStruct hmcResult = null;
          if(useAuxiliaryVariable){
              for (int i = 0; i < nItersPerPathAuxVar; i++)
                  hmcResult = HMC.doIter(rand, L, epsilon, i == 0 ? new DoubleMatrix(initialPoint) : hmcResult.next_q, auxObjective, auxObjective);
          }else{
              hmcResult = HMC.doIter(rand, L, epsilon, new DoubleMatrix(newPoint), noAuxObjective, noAuxObjective);
          }
          newPoint = hmcResult.next_q.data;
      }else
      {
          StopWatch watch = new StopWatch();
          watch.start();
          AHMC ahmc = null;
          if(useAuxiliaryVariable){
              ahmc = AHMC.initializeAHMCWithLBFGS(10000, 1000, auxObjective, auxObjective, initialPoint.length,sizeAdapt);

          }else{
             // we do not recommend without use of auxiliary variable since it can be super slow and cause problems when initializing AHMC
             // ahmc = AHMC.initializeAHMCWithLBFGS(10000, 1000, noAuxObjective, noAuxObjective, initialPoint.length);
              ahmc = new AHMC(10000, 1000, noAuxObjective, noAuxObjective, initialPoint);
          }

          newPoint = ahmc.sample(rand).data;
          epsilon = ahmc.getEpsilon();
          L = ahmc.getL();
          logToFile("time to get epsilon and L in AHMC:" + watch.getTime()/1000);
          logToFile("optimal epsilon" +""+ epsilon);
          logToFile("optimal L" + "" + L);
      }

      mix.parameters.setVector(newPoint);	  
  }
  
  
  private boolean hyperParametersInitialized()
  {
      return epsilon != null;
  }
  
  
  public void logToFile(String someline) {
      this.detailWriter.println(someline);
      this.detailWriter.flush();

  }
  public static ExpectedStatistics<CTMCState> convert(
		  List<PathStatistics> pathStatistics,ExpFamParameters parameters,
		  CTMCStateSpace space)
  {
	  ExpectedStatistics<CTMCState> result = new ExpectedStatistics<CTMCState>(parameters.globalExponentialFamily);
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
}  
