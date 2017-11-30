package conifer.moves;

import java.util.List;
import java.util.Random;

import conifer.RandomUtils.Exponential;
import conifer.RandomUtils.Gamma;
import blang.mcmc.internals.Callback;
import blang.mcmc.MHSampler;
import briefj.collections.UnorderedPair;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.factors.NonClockTreePriorUtils;
import conifer.factors.UnrootedTreeLikelihood;
//import conifer.factors.UnrootedTreeLikelihoodUtils;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.models.MultiCategorySubstitutionModel.PoissonAuxiliarySample;
import conifer.factors.UnrootedTreeLikelihood;



public class AllBranchesScaling extends MHSampler<UnrootedTree>
{
  // original code AllBranchesScaling extends NodeMove
  UnrootedTreeLikelihood likelihood;
  NonClockTreePriorUtils<Exponential.Parameters> prior;
  

  @Override
  public void propose(Random random, Callback callback) {
		// TODO: @Sohrab What does this hack affect?
	    if (random.nextInt(10) != 0)
	      return;
	    
	    // case likelihood to be MultiCategorySubstitutionModel
	    List<PoissonAuxiliarySample> auxVars = ((MultiCategorySubstitutionModel) likelihood.getEvolutionaryModel()).samplePoissonAuxiliaryVariables(random, likelihood.getObservations(), likelihood.getTree());

	    final double alpha = 1.0;
	    final double beta = prior.branchDistributionParameters.getRate();
	    
	    // loop over edges
	    for (UnorderedPair<TreeNode, TreeNode> edge : likelihood.getTree().getTopology().edgeSet())
	    {
	      double updatedShape = alpha;
	      double updatedRate = beta;
	      for (PoissonAuxiliarySample aux : auxVars)
	      {
	        updatedShape += aux.getTransitionCount(edge.getFirst(), edge.getSecond());
	        updatedRate += aux.rate * aux.getSampleCount(edge.getFirst(), edge.getSecond());
	      }
	      // sample new branch length
	      final double newBranchLengths = Gamma.generate(random, updatedRate, updatedShape);
	      variable.updateBranchLength(edge, newBranchLengths);
	    }
  
  } 
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



