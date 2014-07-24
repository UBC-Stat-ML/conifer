package conifer.ctmc.expfam;

import java.io.File;
import java.util.Random;

import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.MCMCRunner;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import briefj.BriefIO;
import briefj.opt.Option;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;


/**
 * Test the phylogenetic MCMC moves on a simple tree model.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class TestRealData extends MCMCRunner
{
  @Option()
  public static RateMtxNames selectedRateMtx = RateMtxNames.POLARITYSIZE;
  
  @DefineFactor(onObservations = true)
  public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = 
    UnrootedTreeLikelihood.fromFastaFile(new File("primates.fasta"), selectedRateMtx).withExpFamMixture(ExpFamMixture.randomGTR(selectedRateMtx));
  
  @DefineFactor
  public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> priorOnParams =
    IIDRealVectorGenerativeFactor.iidNormalOn(likelihood.evolutionaryModel.rateMatrixMixture.parameters);
  
  @DefineFactor
  NonClockTreePrior<RateParameterization> treePrior = NonClockTreePrior.on(likelihood.tree);

  public static void main(String [] args)
  {
    BriefIO.ensureUSLocale();
    TestRealData test = new TestRealData();
    test.factory.mcmcOptions.random = new Random(1243);
    test.run();
  }
}
