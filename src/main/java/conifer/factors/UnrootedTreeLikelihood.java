package conifer.factors;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import bayonet.distributions.Exponential;
import bayonet.marginal.FactorGraph;
import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;
import blang.factors.GenerativeFactor;
import blang.variables.RealVariable;
import briefj.BriefCollections;
import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.ctmc.CTMCParameters;
import conifer.ctmc.RateMatrices;
import conifer.ctmc.SimpleRateMatrix;
import conifer.ctmc.cnv.CopyNumberEmissionModel;
import conifer.ctmc.cnv.CopyNumberMatrix;
import conifer.ctmc.cnv.CopyNumberMixture;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.io.CNParser;
import conifer.io.CopyNumberTreeObservation;
import conifer.io.FastaUtils;
import conifer.io.FixedTreeObservations;
import conifer.io.PhylogeneticObservationFactory;
import conifer.io.TreeObservations;
import conifer.models.CNSpecies;
import conifer.models.DiscreteGammaMixture;
import conifer.models.EvolutionaryModel;
import conifer.models.EvolutionaryModelUtils;
import conifer.models.LikelihoodComputationContext;
import conifer.models.MultiCategorySubstitutionModel;


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
public class UnrootedTreeLikelihood
  <M extends EvolutionaryModel> 
  implements GenerativeFactor
{
  /**
   * 
   */
  @FactorArgument
  public final UnrootedTree tree;
  
  /**
   * 
   */
  @FactorComponent
  public final M evolutionaryModel;
  
  /**
   * 
   */
  @FactorArgument(makeStochastic = true)
  public final TreeObservations observations;
  
  public UnrootedTreeLikelihood(
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
  public static UnrootedTreeLikelihood<MultiCategorySubstitutionModel<DiscreteGammaMixture>> createEmpty(int nSites, List<TreeNode> leaves)
  {
    UnrootedTree tree = defaultTree(leaves);
    SimpleRateMatrix baseRateMatrix = RateMatrices.kimura1980();
    DiscreteGammaMixture gammaMixture = new DiscreteGammaMixture(RealVariable.real(0.1), RealVariable.real(1.0), baseRateMatrix, 4);
    MultiCategorySubstitutionModel<DiscreteGammaMixture> subModel = new MultiCategorySubstitutionModel<DiscreteGammaMixture>(gammaMixture, nSites);
    return new UnrootedTreeLikelihood<MultiCategorySubstitutionModel<DiscreteGammaMixture>>(tree, subModel, new FixedTreeObservations(nSites));
  }
  
  /** 
   * Create a model from a CN file. 
   * @throws IOException 
   * @throws NumberFormatException 
   * 
   * 
   */
  public static UnrootedTreeLikelihood<MultiCategorySubstitutionModel<CopyNumberMixture>> fromCNFile(File cnFile) 
  {
	// read in the raw count data
	Set<CNSpecies> data = CNParser.readCNPairs(cnFile);
        
    // create CNMatrix and CNMixture (gammaMixture not supported)
	CopyNumberMatrix cnMatrix = CopyNumberMatrix.defaultMatrix();
    CopyNumberMixture cnMixture = new CopyNumberMixture(cnMatrix);

    // create treeObservations and initialize the CTMC states
    TreeObservations treeObservations = new CopyNumberTreeObservation(data);
    
    MultiCategorySubstitutionModel<CopyNumberMixture> subModel 
    = new MultiCategorySubstitutionModel<CopyNumberMixture>(cnMixture, treeObservations.nSites());
    
    // make the tree
    UnrootedTree tree = defaultTree(treeObservations.getObservedTreeNodes());
    tree.addNode(TreeNode.withLabel("root"));

    return new UnrootedTreeLikelihood<MultiCategorySubstitutionModel<CopyNumberMixture>>(tree, subModel, treeObservations); 
    
    /*
    // testing (what?)
    Map<String,String> a2s = PhylogeneticObservationFactory.copyNumberFactory().getIndicator2ChunkMap();
 
	StringBuilder result = new StringBuilder();
	for (TreeNode node : observations.getObservedTreeNodes()) {
		result.append(">" + node + "\n");
		double[][] s = (double[][]) observations.get(node);
		String charAtSite = "S";
		for (int j = 0; j < s.length; j++) {
			charAtSite = a2s.get(Arrays.toString(s[j]));
			result.append(charAtSite);
		}
		result.append("\n");
	}
	
	System.out.println(result.toString());
       */ 
  
  }
  
  
  /**
   * Create a model from a FASTA file. Currently assumes nucleotides data and uses kimura 1980 as a default matrix, and a
   * discrete gamma mixture of 4 rates.
   * 
   * @param fastaFile
   * @return
   */
  public static UnrootedTreeLikelihood<MultiCategorySubstitutionModel<DiscreteGammaMixture>> fromFastaFile(File fastaFile)
  {
    Map<TreeNode,CharSequence> data = FastaUtils.readFasta(fastaFile);
    UnrootedTree tree = defaultTree(data.keySet());
    SimpleRateMatrix baseRateMatrix = RateMatrices.kimura1980();
    DiscreteGammaMixture gammaMixture = new DiscreteGammaMixture(RealVariable.real(0.1), RealVariable.real(1.0), baseRateMatrix, 4);
    PhylogeneticObservationFactory factory = PhylogeneticObservationFactory.nucleotidesFactory();
    TreeObservations observations = new FixedTreeObservations(BriefCollections.pick(data.values()).length());
    MultiCategorySubstitutionModel.loadObservations(observations, data, factory);
    
    MultiCategorySubstitutionModel<DiscreteGammaMixture> subModel = new MultiCategorySubstitutionModel<DiscreteGammaMixture>(gammaMixture, observations.nSites());
    return new UnrootedTreeLikelihood<MultiCategorySubstitutionModel<DiscreteGammaMixture>>(tree, subModel, observations);
  }
  
  public static UnrootedTreeLikelihood<MultiCategorySubstitutionModel<DiscreteGammaMixture>> fromFastaProteinFile(File fastaFile)
  {
    Map<TreeNode,CharSequence> data = FastaUtils.readFasta(fastaFile);
    UnrootedTree tree = defaultTree(data.keySet());
    SimpleRateMatrix baseRateMatrix = RateMatrices.accordance();
    DiscreteGammaMixture gammaMixture = new DiscreteGammaMixture(RealVariable.real(0), RealVariable.real(1.0), baseRateMatrix, 1);
    PhylogeneticObservationFactory factory = PhylogeneticObservationFactory.proteinFactory();
    TreeObservations observations = new FixedTreeObservations(BriefCollections.pick(data.values()).length());
    MultiCategorySubstitutionModel.loadObservations(observations, data, factory); 
    MultiCategorySubstitutionModel<DiscreteGammaMixture> subModel = new MultiCategorySubstitutionModel<DiscreteGammaMixture>(gammaMixture, observations.nSites());
    return new UnrootedTreeLikelihood<MultiCategorySubstitutionModel<DiscreteGammaMixture>>(tree, subModel, observations);
  }
  
  public static UnrootedTreeLikelihood<MultiCategorySubstitutionModel<DiscreteGammaMixture>> fromFastaProteinPairFile(File fastaFile)
  {
    Random rand = new Random();
    Map<TreeNode,CharSequence> data = FastaUtils.readFasta(fastaFile);
    UnrootedTree tree = defaultTree(data.keySet());
    SimpleRateMatrix baseRateMatrix = RateMatrices.randomGTR(rand, 400);
    DiscreteGammaMixture gammaMixture = new DiscreteGammaMixture(RealVariable.real(0), RealVariable.real(1.0), baseRateMatrix, 1);
    PhylogeneticObservationFactory factory = PhylogeneticObservationFactory.proteinPairFactory();
    TreeObservations observations = new FixedTreeObservations(factory.nSites()); //;
    //System.out.println(BriefCollections.pick(data.values()).length()/factory.getChunkLength());
    MultiCategorySubstitutionModel.loadObservations(observations, data, factory); 
    MultiCategorySubstitutionModel<DiscreteGammaMixture> subModel = new MultiCategorySubstitutionModel<DiscreteGammaMixture>(gammaMixture, factory.nSites()); //observations.nSites()/factory.getChunkLength());
    return new UnrootedTreeLikelihood<MultiCategorySubstitutionModel<DiscreteGammaMixture>>(tree, subModel, observations);
  }
  
  
  /**
   * Chained method to change the evolutionary model into an exponential family mixture, keeping the other aspects (tree and
   * observations) unchanged.
   * 
   * @param mixture
   * @return
   */
  public UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> withExpFamMixture(ExpFamMixture mixture)
  {
    MultiCategorySubstitutionModel<ExpFamMixture> subModel = new MultiCategorySubstitutionModel<ExpFamMixture>(mixture, this.observations.nSites());
    return new UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>>(this.tree, subModel, this.observations);
  }
  
  public UnrootedTreeLikelihood<MultiCategorySubstitutionModel<DiscreteGammaMixture>> withSingleRateMatrix(CTMCParameters ctmc)
  {
    DiscreteGammaMixture mix = new DiscreteGammaMixture(RealVariable.real(0), RealVariable.real(1.0), ctmc, 1);
    MultiCategorySubstitutionModel<DiscreteGammaMixture> subModel = new MultiCategorySubstitutionModel<DiscreteGammaMixture>(mix, this.observations.nSites()); //observations.nSites()/factory.getChunkLength());
    return new UnrootedTreeLikelihood<MultiCategorySubstitutionModel<DiscreteGammaMixture>>(this.tree, subModel, this.observations);
  }
  
  public UnrootedTreeLikelihood<MultiCategorySubstitutionModel<DiscreteGammaMixture>> withSingleRateMatrix(double [][] matrix)
  {
    return withSingleRateMatrix(new SimpleRateMatrix(matrix, null));
  }
  
  public UnrootedTreeLikelihood<M> withTree(File newickFile)
  {
    return new UnrootedTreeLikelihood<M>(UnrootedTree.fromNewick(newickFile), this.evolutionaryModel, this.observations);
  }
  
  /**
   * 
   * @param leaves
   * @return A default tree obtained by sampling a nonclock tree with an exponential rate one
   *         and a fixed seed.
   */
  public static UnrootedTree defaultTree(Collection<TreeNode> leaves)
  {
    Random rand = new Random(1);
    return NonClockTreePrior.generate(rand, Exponential.on(RealVariable.real()), leaves);
  }

  /**
   * 
   * @return An arbitrary rooting. Does not matter as long as the model is reversible.
   */
  public TreeNode arbitraryNode()
  {
    return TopologyUtils.arbitraryNode(tree);
  }
  
  
  public List<FactorGraph<TreeNode>> buildFactorGraphs()
  {
    return buildFactorGraphs(arbitraryNode());
  }
  
  public List<FactorGraph<TreeNode>> buildFactorGraphs(TreeNode arbitraryRoot)
  {
    return EvolutionaryModelUtils.buildFactorGraphs(evolutionaryModel, tree, arbitraryRoot, observations);
  }

  @Override
  public double logDensity()
  {
    TreeNode arbitraryRoot = arbitraryNode();
    LikelihoodComputationContext context = new LikelihoodComputationContext(buildFactorGraphs(arbitraryRoot), arbitraryRoot);
    return evolutionaryModel.computeLogLikelihood(context);
  }

  @Override
  public void generate(Random random)
  {
    observations.clear();
    if (!observations.getObservedTreeNodes().isEmpty())
      throw new RuntimeException("The method clear() seems to be incorrectly implemented in " + observations.getClass().getName());
    evolutionaryModel.generateObservationsInPlace(random, observations, tree, arbitraryNode());
  }
  
  public static void main(String [] args)
  {
	File f = new File("src/main/resources/conifer/sampleInput/testPatientData.txt");
	UnrootedTreeLikelihood<MultiCategorySubstitutionModel<CopyNumberMixture>> treeLikelihood =
	UnrootedTreeLikelihood.fromCNFile(f);	
	System.out.println(treeLikelihood.observations.toString());
	
	System.err.println("Finished reading CN!");
	  
	  
//    UnrootedTreeLikelihood<MultiCategorySubstitutionModel<DiscreteGammaMixture>> ll = fromFastaFile(new File("/Users/bouchard/Documents/data/utcs/23S.E/R0/cleaned.alignment.fasta"));
//    System.out.println("nSites=" + ll.observations.nSites());
//    System.out.println("nNodes=" + ll.tree.getTopology().vertexSet().size());
//    Random rand = new Random(1);
//    ll.evolutionaryModel.samplePosteriorPaths(rand, ll.observations, ll.tree, BriefCollections.pick(ll.tree.getTopology().vertexSet()), null);
//    System.out.println("done!");
//    ll.
  }
}
