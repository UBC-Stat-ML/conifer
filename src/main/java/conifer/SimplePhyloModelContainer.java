package conifer;

import conifer.ctmc.expfam.ExpFamMixture;
import conifer.factors.UnrootedTreeLikelihoodUtils;
import conifer.models.MultiCategorySubstitutionModel;

public interface SimplePhyloModelContainer {
	
	public UnrootedTreeLikelihoodUtils<MultiCategorySubstitutionModel<ExpFamMixture>> getLikelihood();
}
