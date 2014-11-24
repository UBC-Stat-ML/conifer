package conifer.models;

import org.ejml.simple.SimpleMatrix;

import blang.annotations.FactorArgument;
import conifer.Parsimony;
import conifer.TreeNode;
import conifer.ctmc.CTMC;

/**
 * 
 * @author jewellsean
 *
 * @param <T>
 */

@Deprecated
public class CNMultiCategorySubstitutionModel<T extends RateMatrixMixture> extends
        MultiCategorySubstitutionModel<T>
{

    public CNMultiCategorySubstitutionModel(T rateMatrixMixture, int nSites)
    {
        super(rateMatrixMixture, nSites); 
    }

    // Specify normal rooting
    @Override
    public void buildInitialDistribution(TreeNode node,
        LikelihoodFactoryContext context)
    {
      if (node.toString() != "root")
          throw new RuntimeException("Initial distribution must be built on the root node!");
        
      final int categoryIndex = context.getFactorGraphIndex();
      double[][] ctmc = rateMatrixMixture.getRateMatrix(categoryIndex).getRateMatrix();
      int nStates = ctmc[0].length; 
      double [] stationary = new double[nStates];
      stationary[2] = 1; 
      double [][] allStatios = new double[nSites][stationary.length];
      for (int siteIndex = 0; siteIndex < nSites; siteIndex++)
        allStatios[siteIndex] = stationary;
      context.getDiscreteFactorGraph().unaryTimesEqual(node, new SimpleMatrix(allStatios));
    }
    
    
    
}
