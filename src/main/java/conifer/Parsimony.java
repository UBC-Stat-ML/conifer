package conifer;

import blang.annotations.Processors;
import blang.annotations.Samplers;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.io.CopyNumberTreeObservation;
import conifer.models.CNMultiCategorySubstitutionModel;
import conifer.models.RateMatrixMixture;
import conifer.moves.ParsimonySampler;
import conifer.processors.ParsimonyProcessor;

@Samplers({ParsimonySampler.class})
@Processors({ParsimonyProcessor.class})

public class Parsimony
{

/*
 * Store the min-value for each site in a vector
 * this is the current state of affairs 
 * 
 */
    private ParsimonyVector M; 
    private UnrootedTreeLikelihood<CNMultiCategorySubstitutionModel<RateMatrixMixture>> fullTree;  
    
    public Parsimony(ParsimonyVector M, 
            UnrootedTreeLikelihood<CNMultiCategorySubstitutionModel<RateMatrixMixture>> fullTree)
    {
        this.M = M; 
        this.fullTree = fullTree; 
    }
    
    public UnrootedTreeLikelihood<CNMultiCategorySubstitutionModel<RateMatrixMixture>> getFullTree()
    {
        return fullTree; 
    }
    
    public ParsimonyVector getM()
    {
        return M; 
    }

}
