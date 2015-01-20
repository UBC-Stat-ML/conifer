package conifer.models;

import java.util.Random;

import bayonet.graphs.GraphUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.io.CopyNumberTreeObservation;
import conifer.io.TreeObservations;

/**
 * Neeed to override generative factors 
 * @author jewellsean
 *
 * @param <T>
 */


public class CNMultiCategorySubstitutionModel<T extends RateMatrixMixture> extends
        MultiCategorySubstitutionModel<T>
{

    public CNMultiCategorySubstitutionModel(T rateMatrixMixture, int nSites)
    {
        super(rateMatrixMixture, nSites); 
    }


    @Override
    public void generateObservationsInPlace(
        Random rand, 
        TreeObservations destination,
        UnrootedTree tree, TreeNode root)
    {
      MultiCategoryInternalNodeSample reconstructions = samplePriorInternalNodes(rand, tree, root);
      for (TreeNode leaf : GraphUtils.leaves(tree.getTopology()))
        {
          ((CopyNumberTreeObservation) destination).set(rand, leaf, reconstructions.internalIndicators.get(leaf));
        }
    }
    
    
    
}
