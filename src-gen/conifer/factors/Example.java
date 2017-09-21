package conifer.factors;

import blang.core.ConstantSupplier;
import blang.core.DeboxedName;
import blang.core.Model;
import blang.core.ModelBuilder;
import blang.core.ModelComponent;
import blang.core.Param;
import blang.inits.Arg;
import blang.types.Index;
import blang.types.Plate;
import blang.types.Plated;
import ca.ubc.stat.blang.StaticJavaUtils;
import conifer.UnrootedTree;
import conifer.factors.RealValuedDensity;
import conifer.io.TreeObservations;
import conifer.models.EvolutionaryModel;
import java.util.ArrayList;
import java.util.Collection;
import java.util.function.Supplier;

@SuppressWarnings("all")
public class Example implements Model {
  public static class Builder implements ModelBuilder {
    @Arg
    public Plate<String> datasets;
    
    @Arg
    public Plated<TreeObservations> observations;
    
    @Arg
    public Plated<UnrootedTree> trees;
    
    @Arg
    public EvolutionaryModel evolutionaryModel;
    
    @Arg
    public RealValuedDensity branchLengthDensity;
    
    public Example build() {
      // For each optional type, either get the value, or evaluate the ?: expression
      final Plate<String> __datasets = datasets;
      final Plated<TreeObservations> __observations = observations;
      final Plated<UnrootedTree> __trees = trees;
      final EvolutionaryModel __evolutionaryModel = evolutionaryModel;
      final RealValuedDensity __branchLengthDensity = branchLengthDensity;
      // Build the instance after boxing params
      return new Example(
        __observations, 
        __trees, 
        new ConstantSupplier(__datasets), 
        new ConstantSupplier(__evolutionaryModel), 
        new ConstantSupplier(__branchLengthDensity)
      );
    }
  }
  
  @Param
  private final Supplier<Plate<String>> $generated__datasets;
  
  public Plate<String> getDatasets() {
    return $generated__datasets.get();
  }
  
  private final Plated<TreeObservations> observations;
  
  public Plated<TreeObservations> getObservations() {
    return observations;
  }
  
  private final Plated<UnrootedTree> trees;
  
  public Plated<UnrootedTree> getTrees() {
    return trees;
  }
  
  @Param
  private final Supplier<EvolutionaryModel> $generated__evolutionaryModel;
  
  public EvolutionaryModel getEvolutionaryModel() {
    return $generated__evolutionaryModel.get();
  }
  
  @Param
  private final Supplier<RealValuedDensity> $generated__branchLengthDensity;
  
  public RealValuedDensity getBranchLengthDensity() {
    return $generated__branchLengthDensity.get();
  }
  
  /**
   * Utility main method for posterior inference on this model
   */
  public static void main(final String[] arguments) {
    StaticJavaUtils.callRunner(Builder.class, arguments);
  }
  
  /**
   * Auxiliary method generated to translate:
   * datasets.indices
   */
  private static Iterable<Index<String>> $generated__0(final Plate<String> datasets, final Plated<TreeObservations> observations, final Plated<UnrootedTree> trees, final EvolutionaryModel evolutionaryModel, final RealValuedDensity branchLengthDensity) {
    Iterable<Index<String>> _indices = datasets.indices();
    return _indices;
  }
  
  /**
   * Auxiliary method generated to translate:
   * trees.get(dataset)
   */
  private static UnrootedTree $generated__1(final Index<String> dataset, final Plate<String> datasets, final Plated<TreeObservations> observations, final Plated<UnrootedTree> trees, final EvolutionaryModel evolutionaryModel, final RealValuedDensity branchLengthDensity) {
    UnrootedTree _get = trees.get(dataset);
    return _get;
  }
  
  /**
   * Auxiliary method generated to translate:
   * branchLengthDensity
   */
  private static RealValuedDensity $generated__2(final RealValuedDensity branchLengthDensity) {
    return branchLengthDensity;
  }
  
  public static class $generated__2_class implements Supplier<RealValuedDensity> {
    public RealValuedDensity get() {
      return $generated__2($generated__branchLengthDensity.get());
    }
    
    public String toString() {
      return "branchLengthDensity";
    }
    
    private Supplier<RealValuedDensity> $generated__branchLengthDensity;
    
    public $generated__2_class(final Supplier<RealValuedDensity> $generated__branchLengthDensity) {
      this.$generated__branchLengthDensity = $generated__branchLengthDensity;
    }
  }
  
  /**
   * Auxiliary method generated to translate:
   * trees.get(dataset)
   */
  private static UnrootedTree $generated__3(final Index<String> dataset, final Plate<String> datasets, final Plated<TreeObservations> observations, final Plated<UnrootedTree> trees, final EvolutionaryModel evolutionaryModel, final RealValuedDensity branchLengthDensity) {
    UnrootedTree _get = trees.get(dataset);
    return _get;
  }
  
  /**
   * Auxiliary method generated to translate:
   * observations.get(dataset)
   */
  private static TreeObservations $generated__4(final Index<String> dataset, final Plate<String> datasets, final Plated<TreeObservations> observations, final Plated<UnrootedTree> trees, final EvolutionaryModel evolutionaryModel, final RealValuedDensity branchLengthDensity) {
    TreeObservations _get = observations.get(dataset);
    return _get;
  }
  
  /**
   * Auxiliary method generated to translate:
   * tree
   */
  private static UnrootedTree $generated__5(final UnrootedTree tree, final EvolutionaryModel evolutionaryModel) {
    return tree;
  }
  
  public static class $generated__5_class implements Supplier<UnrootedTree> {
    public UnrootedTree get() {
      return $generated__5(tree, $generated__evolutionaryModel.get());
    }
    
    public String toString() {
      return "tree";
    }
    
    private UnrootedTree tree;
    
    private Supplier<EvolutionaryModel> $generated__evolutionaryModel;
    
    public $generated__5_class(final UnrootedTree tree, final Supplier<EvolutionaryModel> $generated__evolutionaryModel) {
      this.tree = tree;
      this.$generated__evolutionaryModel = $generated__evolutionaryModel;
    }
  }
  
  /**
   * Auxiliary method generated to translate:
   * evolutionaryModel
   */
  private static EvolutionaryModel $generated__6(final UnrootedTree tree, final EvolutionaryModel evolutionaryModel) {
    return evolutionaryModel;
  }
  
  public static class $generated__6_class implements Supplier<EvolutionaryModel> {
    public EvolutionaryModel get() {
      return $generated__6(tree, $generated__evolutionaryModel.get());
    }
    
    public String toString() {
      return "evolutionaryModel";
    }
    
    private UnrootedTree tree;
    
    private Supplier<EvolutionaryModel> $generated__evolutionaryModel;
    
    public $generated__6_class(final UnrootedTree tree, final Supplier<EvolutionaryModel> $generated__evolutionaryModel) {
      this.tree = tree;
      this.$generated__evolutionaryModel = $generated__evolutionaryModel;
    }
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
  public Example(@DeboxedName("observations") final Plated<TreeObservations> observations, @DeboxedName("trees") final Plated<UnrootedTree> trees, @Param @DeboxedName("datasets") final Supplier<Plate<String>> $generated__datasets, @Param @DeboxedName("evolutionaryModel") final Supplier<EvolutionaryModel> $generated__evolutionaryModel, @Param @DeboxedName("branchLengthDensity") final Supplier<RealValuedDensity> $generated__branchLengthDensity) {
    this.$generated__datasets = $generated__datasets;
    this.observations = observations;
    this.trees = trees;
    this.$generated__evolutionaryModel = $generated__evolutionaryModel;
    this.$generated__branchLengthDensity = $generated__branchLengthDensity;
  }
  
  /**
   * A component can be either a distribution, support constraint, or another model  
   * which recursively defines additional components.
   */
  public Collection<ModelComponent> components() {
    ArrayList<ModelComponent> components = new ArrayList();
    
    for (Index<String> dataset : $generated__0($generated__datasets.get(), observations, trees, $generated__evolutionaryModel.get(), $generated__branchLengthDensity.get())) {
      { // Code generated by: trees.get(dataset) | branchLengthDensity ~ NonClockTreePrior(branchLengthDensity)
        // Construction and addition of the factor/model:
        components.add(
          new conifer.factors.NonClockTreePrior(
            $generated__1(dataset, $generated__datasets.get(), observations, trees, $generated__evolutionaryModel.get(), $generated__branchLengthDensity.get()), 
            new $generated__2_class($generated__branchLengthDensity)
          )
          );
      }
      { // Code generated by: observations.get(dataset) | UnrootedTree tree = trees.get(dataset), evolutionaryModel ~ UnrootedTreeLikelihood(tree, evolutionaryModel)
        // Required initialization:
        UnrootedTree tree = $generated__3(dataset, $generated__datasets.get(), observations, trees, $generated__evolutionaryModel.get(), $generated__branchLengthDensity.get());
        // Construction and addition of the factor/model:
        components.add(
          new conifer.factors.UnrootedTreeLikelihood(
            $generated__4(dataset, $generated__datasets.get(), observations, trees, $generated__evolutionaryModel.get(), $generated__branchLengthDensity.get()), 
            new $generated__5_class(tree, $generated__evolutionaryModel), 
            new $generated__6_class(tree, $generated__evolutionaryModel)
          )
          );
      }
    }
    
    return components;
  }
}
