package conifer;

import conifer.ctmc.expfam.ExpFamMixture;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;

public interface SimplePhyloModelContainer {
	
	public UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> getLikelihood();
}
