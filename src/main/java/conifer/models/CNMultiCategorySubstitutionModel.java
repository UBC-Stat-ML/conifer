package conifer.models;

import blang.annotations.FactorComponent;

/**
 * 
 * @author jewellsean
 *
 * @param <T>
 */
public class CNMultiCategorySubstitutionModel<T extends RateMatrixMixture> extends
        MultiCategorySubstitutionModel<T>
{

    public CNMultiCategorySubstitutionModel(T rateMatrixMixture, ParsimonyModel pm, int nSites)
    {
        super(rateMatrixMixture, nSites);
        this.pm = pm; 
    }

    @FactorComponent
    public final ParsimonyModel pm; 
    
}
