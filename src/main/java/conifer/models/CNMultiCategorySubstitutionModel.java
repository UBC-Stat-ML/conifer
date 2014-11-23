package conifer.models;

import blang.annotations.FactorArgument;
import conifer.Parsimony;

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

    public CNMultiCategorySubstitutionModel(T rateMatrixMixture, Parsimony parsimony, int nSites)
    {
        super(rateMatrixMixture, nSites);
        this.parsimony = parsimony; 
    }

    public final Parsimony parsimony; 
    
}
