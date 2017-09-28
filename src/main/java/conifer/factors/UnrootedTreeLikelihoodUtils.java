package conifer.factors;

import java.io.File;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Random;

import conifer.RandomUtils.Exponential;
import bayonet.marginal.FactorGraph;
import blang.types.RealScalar;
import briefj.BriefCollections;
import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.ctmc.CTMCParameters;
import conifer.ctmc.RateMatrices;
import conifer.ctmc.SimpleRateMatrix;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.io.FastaUtils;
import conifer.io.FixedTreeObservations;
import conifer.io.PhylogeneticObservationFactory;
import conifer.io.TreeObservations;
import conifer.models.DiscreteGammaMixture;
import conifer.models.EvolutionaryModel;
import conifer.models.EvolutionaryModelUtils;
import conifer.models.LikelihoodComputationContext;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.ctmc.expfam.RateMtxNames;

/**
 * Many phylogenetic likelihoods can be viewed as a factor that connects a 
 * tree, an evolutionary model, and observations.
 * 
 * In this context, an evolutionary model is basically a factory that builds 
 * a factor graph given a tree and observations (actually, a list of factor 
 * graphs to handle multicategory/PIP/multiple genes).
 * 
 * Note that the main assumption we are making here is that the model is assumed to 
 * be reversible (we do not keep track of a root).
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 * @param <M>
 */
public class UnrootedTreeLikelihoodUtils<M>
{
  /**
   * 
   */
  // TODO: @Sohrab, is this the best way to supress sampling of the tree (fixed topology and branch length?)
  public final UnrootedTree tree;
  
  /**
   * 
   */

  public final M evolutionaryModel;
  
  /**
   * 
   */
  public final TreeObservations observations;
  
  public UnrootedTreeLikelihoodUtils(
      UnrootedTree tree, 
      M evolutionaryModel,
      TreeObservations observations)
  {
    this.tree = tree;
    this.evolutionaryModel = evolutionaryModel;
    this.observations = observations;
  }
  
  /**
   * Creates a model with nSites sites, default values, and 
   * no data. Used mostly for simulation-based testing.
   * @return
   */
  public static  UnrootedTreeLikelihoodUtils<MultiCategorySubstitutionModel<DiscreteGammaMixture>> createEmpty(int nSites, List<TreeNode> leaves, RateMtxNames selectedRateMtx)
  {
    UnrootedTree tree = defaultTree(leaves);
    SimpleRateMatrix baseRateMatrix = RateMatrices.rateMtxModel(selectedRateMtx);
    DiscreteGammaMixture gammaMixture = new DiscreteGammaMixture(new RealScalar(0.1), new RealScalar(1.0), baseRateMatrix, 4);
    MultiCategorySubstitutionModel<DiscreteGammaMixture> subModel = new MultiCategorySubstitutionModel<DiscreteGammaMixture>(gammaMixture, nSites);
    return new UnrootedTreeLikelihoodUtils<MultiCategorySubstitutionModel<DiscreteGammaMixture>>(tree, subModel, new FixedTreeObservations(nSites));
  }

  /**
   * Creates a model with nSites sites, default values and fixed tree topology from reading a newick tree file
   * @param nSites
   * @param treeFile
   * @param selectedRateMtx
   * @return
   */

  public static  UnrootedTreeLikelihoodUtils<MultiCategorySubstitutionModel<DiscreteGammaMixture>> createEmptyWithFixedTree(int nSites, File treeFile, RateMtxNames selectedRateMtx)
  {
    UnrootedTree tree = UnrootedTree.fromNewick(treeFile);
    SimpleRateMatrix baseRateMatrix = RateMatrices.rateMtxModel(selectedRateMtx);
    DiscreteGammaMixture gammaMixture = new DiscreteGammaMixture(new RealScalar(0), new RealScalar(1.0), baseRateMatrix, 1);
    MultiCategorySubstitutionModel<DiscreteGammaMixture> subModel = new MultiCategorySubstitutionModel<DiscreteGammaMixture>(gammaMixture, nSites);
    return new UnrootedTreeLikelihoodUtils<MultiCategorySubstitutionModel<DiscreteGammaMixture>>(tree, subModel, new FixedTreeObservations(nSites));
  }
  
  /**
   * Create a model from a FASTA file. Currently assumes nucleotides data and uses kimura 1980 as a default matrix, and a
   * discrete gamma mixture of 4 rates.
   * 
   * @param fastaFile
   * @return
   */
  public static UnrootedTreeLikelihoodUtils<MultiCategorySubstitutionModel<DiscreteGammaMixture>> fromFastaFile(File fastaFile, final RateMtxNames selectedRateMtx)
  {
    Map<TreeNode,CharSequence> data = FastaUtils.readFasta(fastaFile);
    UnrootedTree tree = defaultTree(data.keySet());
    SimpleRateMatrix baseRateMatrix = RateMatrices.rateMtxModel(selectedRateMtx);
    DiscreteGammaMixture gammaMixture = new DiscreteGammaMixture(new RealScalar(0.1),new RealScalar(1.0), baseRateMatrix, 4);
    PhylogeneticObservationFactory factory = PhylogeneticObservationFactory.selectedFactory(selectedRateMtx);
    TreeObservations observations = new FixedTreeObservations(BriefCollections.pick(data.values()).length()/factory.getChunkLength());
    MultiCategorySubstitutionModel.loadObservations(observations, data, factory); 
    MultiCategorySubstitutionModel<DiscreteGammaMixture> subModel = new MultiCategorySubstitutionModel<DiscreteGammaMixture>(gammaMixture, observations.nSites());
    return new UnrootedTreeLikelihoodUtils<MultiCategorySubstitutionModel<DiscreteGammaMixture>>(tree, subModel, observations);
  }
  
    
 
//  
//  /**
//   * 
//   * @param leaves
//   * @return A default tree obtained by sampling a nonclock tree with an exponential rate one
//   *         and a fixed seed.
//   */
//  public static UnrootedTree defaultTree(Collection<TreeNode> leaves)
//  {
//    Random rand = new Random(1);
//    return NonClockTreePriorUtils.sample(rand, Exponential.on(new RealScalar(0)), leaves);
//  }

  
  public static void main(String [] args)
  {
    RateMtxNames selectedRateMtx = RateMtxNames.KIMURA1980;
    UnrootedTreeLikelihoodUtils<MultiCategorySubstitutionModel<DiscreteGammaMixture>> ll = fromFastaFile(new File("/Users/bouchard/Documents/data/utcs/23S.E/R0/cleaned.alignment.fasta"), selectedRateMtx);
    System.out.println("nSites=" + ll.observations.nSites());
    System.out.println("nNodes=" + ll.tree.getTopology().vertexSet().size());
    Random rand = new Random(1);
    ll.evolutionaryModel.samplePosteriorPaths(rand, ll.observations, ll.tree, BriefCollections.pick(ll.tree.getTopology().vertexSet()), null);
    System.out.println("done!");
//    
  }
}
