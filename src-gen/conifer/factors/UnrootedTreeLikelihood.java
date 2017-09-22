package conifer.factors;

import bayonet.marginal.FactorGraph;
import blang.core.ConstantSupplier;
import blang.core.DeboxedName;
import blang.core.ForwardSimulator;
import blang.core.LogScaleFactor;
import blang.core.Model;
import blang.core.ModelBuilder;
import blang.core.ModelComponent;
import blang.core.Param;
import blang.inits.Arg;
import ca.ubc.stat.blang.StaticJavaUtils;
import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.io.TreeObservations;
import conifer.models.EvolutionaryModel;
import conifer.models.EvolutionaryModelUtils;
import conifer.models.LikelihoodComputationContext;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Random;
import java.util.function.Supplier;

@SuppressWarnings("all")
public class UnrootedTreeLikelihood implements Model, ForwardSimulator {
  public static class Builder implements ModelBuilder {
    @Arg
    public TreeObservations observations;
    
    @Arg
    public UnrootedTree tree;
    
    @Arg
    public EvolutionaryModel evolutionaryModel;
    
    public UnrootedTreeLikelihood build() {
      // For each optional type, either get the value, or evaluate the ?: expression
      final TreeObservations __observations = observations;
      final UnrootedTree __tree = tree;
      final EvolutionaryModel __evolutionaryModel = evolutionaryModel;
      // Build the instance after boxing params
      return new UnrootedTreeLikelihood(
        __observations, 
        new ConstantSupplier(__tree), 
        new ConstantSupplier(__evolutionaryModel)
      );
    }
  }
  
  private final TreeObservations observations;
  
  public TreeObservations getObservations() {
    return observations;
  }
  
  @Param
  private final Supplier<UnrootedTree> $generated__tree;
  
  public UnrootedTree getTree() {
    return $generated__tree.get();
  }
  
  @Param
  private final Supplier<EvolutionaryModel> $generated__evolutionaryModel;
  
  public EvolutionaryModel getEvolutionaryModel() {
    return $generated__evolutionaryModel.get();
  }
  
  /**
   * Utility main method for posterior inference on this model
   */
  public static void main(final String[] arguments) {
    StaticJavaUtils.callRunner(Builder.class, arguments);
  }
  
  /**
   * Auxiliary method generated to translate:
   * { val TreeNode arbitraryRoot = TopologyUtils.arbitraryNode(tree) val LikelihoodComputationContext context = new LikelihoodComputationContext(EvolutionaryModelUtils.buildFactorGraphs(evolutionaryModel, tree, arbitraryRoot, observations), arbitraryRoot) return evolutionaryModel.computeLogLikelihood(context) }
   */
  private static double $generated__0(final UnrootedTree tree, final TreeObservations observations, final EvolutionaryModel evolutionaryModel) {
    final TreeNode arbitraryRoot = TopologyUtils.arbitraryNode(tree);
    List<FactorGraph<TreeNode>> _buildFactorGraphs = EvolutionaryModelUtils.buildFactorGraphs(evolutionaryModel, tree, arbitraryRoot, observations);
    final LikelihoodComputationContext context = new LikelihoodComputationContext(_buildFactorGraphs, arbitraryRoot);
    return evolutionaryModel.computeLogLikelihood(context);
  }
  
  public static class $generated__0_class implements LogScaleFactor {
    public double logDensity() {
      return $generated__0($generated__tree.get(), observations, $generated__evolutionaryModel.get());
    }
    
    public String toString() {
      return "{ val TreeNode arbitraryRoot = TopologyUtils.arbitraryNode(tree) val LikelihoodComputationContext context = new LikelihoodComputationContext(EvolutionaryModelUtils.buildFactorGraphs(evolutionaryModel, tree, arbitraryRoot, observations), arbitraryRoot) return evolutionaryModel.computeLogLikelihood(context) }";
    }
    
    private Supplier<UnrootedTree> $generated__tree;
    
    private TreeObservations observations;
    
    private Supplier<EvolutionaryModel> $generated__evolutionaryModel;
    
    public $generated__0_class(final Supplier<UnrootedTree> $generated__tree, final TreeObservations observations, final Supplier<EvolutionaryModel> $generated__evolutionaryModel) {
      this.$generated__tree = $generated__tree;
      this.observations = observations;
      this.$generated__evolutionaryModel = $generated__evolutionaryModel;
    }
  }
  
  /**
   * Auxiliary method generated to translate:
   * { observations.clear() if (!observations.getObservedTreeNodes().isEmpty()) throw new RuntimeException("The method clear() seems to be incorrectly implemented in " + observations.getClass().getName()) evolutionaryModel.generateObservationsInPlace(rand, observations, tree, TopologyUtils.arbitraryNode(tree)) }
   */
  private static void $generated__1(final Random rand, final TreeObservations observations, final UnrootedTree tree, final EvolutionaryModel evolutionaryModel) {
    observations.clear();
    List<TreeNode> _observedTreeNodes = observations.getObservedTreeNodes();
    boolean _isEmpty = _observedTreeNodes.isEmpty();
    boolean _not = (!_isEmpty);
    if (_not) {
      Class<? extends TreeObservations> _class = observations.getClass();
      String _name = _class.getName();
      String _plus = ("The method clear() seems to be incorrectly implemented in " + _name);
      throw new RuntimeException(_plus);
    }
    TreeNode _arbitraryNode = TopologyUtils.arbitraryNode(tree);
    evolutionaryModel.generateObservationsInPlace(rand, observations, tree, _arbitraryNode);
  }
  
  /**
   * Note: the generated code has the following properties used at runtime:
   *   - all arguments are annotated with a BlangVariable annotation
   *   - params additionally have a Param annotation
   *   - the order of the arguments is as follows:
   *     - first, all the random variables in the order they occur in the blang file
   *     - second, all the params in the order they occur in the blang file
   * 
   */
  public UnrootedTreeLikelihood(@DeboxedName("observations") final TreeObservations observations, @Param @DeboxedName("tree") final Supplier<UnrootedTree> $generated__tree, @Param @DeboxedName("evolutionaryModel") final Supplier<EvolutionaryModel> $generated__evolutionaryModel) {
    this.observations = observations;
    this.$generated__tree = $generated__tree;
    this.$generated__evolutionaryModel = $generated__evolutionaryModel;
  }
  
  /**
   * A component can be either a distribution, support constraint, or another model  
   * which recursively defines additional components.
   */
  public Collection<ModelComponent> components() {
    ArrayList<ModelComponent> components = new ArrayList();
    
    { // Code generated by: (tree, observations, evolutionaryModel) { val TreeNode arbitraryRoot = TopologyUtils.arbitraryNode(tree) val LikelihoodComputationContext context = new LikelihoodComputationContext(EvolutionaryModelUtils.buildFactorGraphs(evolutionaryModel, tree, arbitraryRoot, observations), arbitraryRoot) return evolutionaryModel.computeLogLikelihood(context) }
      // Construction and addition of the factor/model:
      components.add(
        new $generated__0_class($generated__tree, observations, $generated__evolutionaryModel));
    }
    
    return components;
  }
  
  public void generate(final Random rand) {
    $generated__1(rand, observations, $generated__tree.get(), $generated__evolutionaryModel.get());
  }
}
