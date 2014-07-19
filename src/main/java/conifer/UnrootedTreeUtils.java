package conifer;


import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;


import org.apache.commons.lang3.tuple.Pair;
import org.jgrapht.Graphs;
import org.jgrapht.UndirectedGraph;


import bayonet.graphs.GraphUtils;
import bayonet.marginal.BinaryFactor;
import bayonet.marginal.FactorGraph;
import bayonet.marginal.FactorOperations;
import bayonet.marginal.UnaryFactor;
import bayonet.marginal.algo.SumProduct;
import briefj.BriefCollections;
import briefj.BriefIO;
import briefj.BriefLists;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.Tree;
import briefj.collections.UnorderedPair;

import com.beust.jcommander.internal.Maps;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import conifer.io.newick.NewickParser;
import conifer.io.newick.ParseException;



public class UnrootedTreeUtils
{
  /**
   * See comments in UnrootedTree
   * @param f
   * @return
   */
  public static UnrootedTree fromNewick(File f)
  {
    return fromNewickString(BriefIO.fileToString(f));
  }
  
  /**
   * See comments in UnrootedTree
   * @param f
   * @return
   */
  public static UnrootedTree fromNewickString(String newickString)
  {
    NewickParser parser = new NewickParser(newickString); 
    try
    {
      Tree<TreeNode> topo = parser.parse();
      Map<TreeNode,Double> bls = parser.getBranchLengths();
      UnrootedTree result = new UnrootedTree();
      process(result, topo, bls, null);
      return result;
    } catch (ParseException e)
    {
      throw new RuntimeException(e);
    }
  }
  
  private static void process(
      UnrootedTree result, Tree<TreeNode> topo,
      Map<TreeNode, Double> bls, TreeNode ancestor)
  {
    TreeNode current = topo.getLabel();
    result.addNode(current);
    if (ancestor != null)
      result.addEdge(ancestor, current, bls.get(current));
    for (Tree<TreeNode> child : topo.getChildren())
      process(result, child, bls, current);
  }

  /**
   * Creates a newick string for the given tree. 
   * @param tree
   * @return
   */
  public static String toNewick(UnrootedTree tree)
  {
    StringBuilder result = new StringBuilder();
    TreeNode pseudoRoot = null;
    // root the tree at an internal node if possible
    List<TreeNode> internalNodes = GraphUtils.internalNodes(tree.getTopology());
    if (!internalNodes.isEmpty())
      pseudoRoot = BriefCollections.pick(internalNodes);
    else if (!tree.getTopology().vertexSet().isEmpty()) // otherwise, do it arbitrarily
      pseudoRoot = BriefCollections.pick(tree.getTopology().vertexSet());
    if (pseudoRoot != null)
      toNewick(tree, null, pseudoRoot, result);
    result.append(";");
    return result.toString();
  }
  
  private static void toNewick(UnrootedTree tree, TreeNode parent, TreeNode current, StringBuilder builder)
  {
    List<TreeNode> children = Graphs.neighborListOf(tree.getTopology(), current);
    if (parent != null)
      if (!children.remove(parent))
        throw new RuntimeException();
    
    if (!children.isEmpty())
    {
      builder.append("(");
      for (int cIndex = 0; cIndex < children.size(); cIndex++)
      {
        toNewick(tree, current, children.get(cIndex), builder);
        if (cIndex != children.size() - 1)
          builder.append(",");
      }
      builder.append(")");
    }
    if (current.isLabelled())
    {
      String label = current.toString();
      if (label.contains("(") || 
          label.contains(")") || 
          label.contains(",") || 
          label.contains(":") ||
          label.contains(";"))
        throw new RuntimeException();
      builder.append(label);
    }
    if (parent != null)
      builder.append(":" + tree.getBranchLength(current, parent));
  }
  
  /**
   * Compute the total branch length distance between each pair of leaves in a tree.
   * 
   * By total branch length distance, we mean the sum of branch lengths encountered
   * in the unique shortest path between two leaves.
   * 
   * @param tree
   * @return
   */
  public static Counter<UnorderedPair<TreeNode,TreeNode>> allTotalBranchLengthDistances(UnrootedTree tree)
  {
    return new EfficientUnrootedTree(tree).allTotalBranchLengthDistances();
  }
  
  /**
   * An efficient array-based implementation of unrooted trees used
   * internally to compute pairwise distances.
   * 
   * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
   *
   */
  private static class EfficientUnrootedTree
  {
    private final int [][] nbhrs;
    private final double [][] bls;
    private final int nLeaves;
    private final Indexer<TreeNode> indexer = new Indexer<TreeNode>();
    
    public EfficientUnrootedTree(UnrootedTree ut)
    {
      final List<TreeNode> leaves = ut.leaves();
      this.nLeaves = leaves.size();
      for (final TreeNode leaf : leaves)
        indexer.addToIndex(leaf);
      for (final TreeNode node : ut.getTopology().vertexSet())
        if (!indexer.containsObject(node))
          indexer.addToIndex(node);
      final int size = indexer.size();
      this.nbhrs = new int[size][];
      this.bls = new double[size][];
      for (final TreeNode t1 : ut.getTopology().vertexSet())
      {
        final int i1 = indexer.o2i(t1);
        int currentDegree = ut.getTopology().degreeOf(t1);
        nbhrs[i1] = new int[currentDegree];
        bls[i1] = new double[currentDegree];
        int idx = 0;
        for (TreeNode t2 : Graphs.neighborListOf(ut.getTopology(), t1))
        {
          final int i2 = indexer.o2i(t2);
          nbhrs[i1][idx] = i2;
          bls[i1][idx] = ut.getBranchLength(t1, t2);
          idx++;
        }
      }
    }
    
    public Counter<UnorderedPair<TreeNode,TreeNode>> allTotalBranchLengthDistances()
    {
      double [][] result = new double[nLeaves][nLeaves];
      
      for (int start = 0; start < nLeaves - 1; start++)
        _dfsTotalBL(0.0, start, start, -1, result);
      
      // conversion
      Counter<UnorderedPair<TreeNode,TreeNode>> convertedResult = new Counter<UnorderedPair<TreeNode,TreeNode>>();
      for (int l1 = 0; l1 < nLeaves; l1++)
      {
        final TreeNode t1 = indexer.i2o(l1);
        for (int l2 = l1 + 1; l2 < nLeaves; l2++)
          convertedResult.setCount(new UnorderedPair<TreeNode, TreeNode>(t1, indexer.i2o(l2)), result[l1][l2]);
      }
      return convertedResult;
    }

    private void _dfsTotalBL(double parentLen, int start, int current, int parent, double[][] result)
    {
      final int [] thisNbhrs = nbhrs[current];
      final double [] thisBLs = bls[current];
      if (thisNbhrs.length != thisBLs.length)
        throw new RuntimeException();
      for (int nIndex = 0; nIndex < thisNbhrs.length; nIndex++)
      {
        int nbhr = thisNbhrs[nIndex];
        if (nbhr != parent)
        {
          final double newLen = parentLen + thisBLs[nIndex];
          if (nbhr < nLeaves)
            // a leaf!
            result[start][nbhr] = newLen;
          else
            // recurse!
            _dfsTotalBL(newLen, start, nbhr, current, result);
        }
      }
    }
  }
  
  /**
   * Remove nodes that have exactly two neighbors, changing branch lengths appropriately 
   * @param t
   */
  public static void simplify(UnrootedTree t)
  {
    simplify(t, t.leaves().get(0), null);
  }
  
  /**
   * See comments in unrooted
   * 
   * @param t
   * @param parent
   * @param current
   * @param accumulatedBranchLength
   * @return
   */
  private static void simplify(UnrootedTree t, TreeNode goodNode, TreeNode parent)
  {
    outerLoop : for (UnorderedPair<TreeNode,TreeNode> neighborEdge : Lists.newArrayList(t.getTopology().edgesOf(goodNode)))
    {
      TreeNode neighborNode = GraphUtils.pickOther(neighborEdge, goodNode);
      if (neighborNode == parent)
        continue outerLoop;
      
      // follow down until we get to a good node,
      TreeNode previous = goodNode;
      TreeNode current = neighborNode;
      double totalBranchLength = 0.0;
      List<Pair<TreeNode,TreeNode>> path = Lists.newArrayList();
      do
      {
        totalBranchLength += t.getBranchLength(previous, current);
        path.add(Pair.of(previous, current));
        TreeNode next = findNext(t.getTopology().edgesOf(current), current, previous);
        previous = current;
        current = next;
      } while (previous != null && t.getTopology().edgesOf(previous).size() == 2);
      
      // if there were more than one hop to a good node, do simplifications
      TreeNode last = BriefLists.last(path).getRight();
      if (path.size() > 1)
      {
        for (int i = 0; i < path.size(); i++)
        {
          Pair<TreeNode,TreeNode> currentEdge = path.get(i);
          t.removeEdge(currentEdge.getLeft(), currentEdge.getRight());
          if (i != 0)
            t.getTopology().removeVertex(currentEdge.getLeft());
        }
        t.addEdge(goodNode, last, totalBranchLength);
      }
      
      simplify(t, last, goodNode);
    }

  }
  private static TreeNode findNext(
      Set<UnorderedPair<TreeNode, TreeNode>> edges, TreeNode center, TreeNode previous)
  {
    for (UnorderedPair<TreeNode, TreeNode> edge : edges)
    {
      TreeNode node = GraphUtils.pickOther(edge, center);
      if (node != previous)
        return node;
    }
    return null;
  }

  /**
   * The sum of all branch lengths.
   * 
   * @param tree
   * @return
   */
  public static double totalTreeLength(UnrootedTree tree)
  {
    double sum = 0.0;
    for (double branchLength : tree.getBranchLengths().values())
      sum += branchLength;
    return sum;
  }
  
  /**
   * Compute l0 (topology only), l1, and l2 metrics between the two provided trees.
   * 
   * @param truth
   * @param guess
   * @return
   */
  public static Map<TreeMetric,Double> computeTreeMetrics(UnrootedTree truth, UnrootedTree guess)
  {
    Indexer<TreeNode> tipsIndexer = tipIndexer(truth);
    CladeCalculator 
      calc1 = cladeCalculator(truth, tipsIndexer),
      calc2 = cladeCalculator(guess, tipsIndexer);
    Map<TreeMetric,Double> result = new HashMap<TreeMetric, Double>();
    
    double l0 = 0.0, l1 = 0.0, l2 = 0.0;
    Counter<UnorderedPair<Clade,Clade>> 
      clades1 = calc1.bipartitions(),
      clades2 = calc2.bipartitions();
    
    for (UnorderedPair<Clade,Clade> bipart : Sets.union(clades1.keySet(), clades2.keySet()))
    {
      final double 
        b1 = clades1.getCount(bipart),
        b2 = clades2.getCount(bipart);
      l0 += l0(b1,b2);
      l1 += l1(b1,b2);
      l2 += l2(b1,b2);
    }
    
    result.put(TreeMetric.l0, l0);
    result.put(TreeMetric.l1, l1);
    result.put(TreeMetric.l2, Math.sqrt(l2));
    
    return result;
  }
  
  private static double l2(double b1, double b2)
  {
    final double x = b1 - b2;
    return x * x;
  }

  private static double l1(double b1, double b2)
  {
    return Math.abs(b1 - b2);
  }

  private static double l0(double b1, double b2)
  {
    return Math.abs(isPositive(b1) - isPositive(b2));
  }

  private static double isPositive(double b2)
  {
    return b2 > 0 ? 1.0 : 0.0;
  }

  public static enum TreeMetric { l0, l1, l2 };
  
  
  public static double bipartitionTopologySymmetricDifference(CladeCalculator calc1, CladeCalculator calc2)
  {
    Set<UnorderedPair<Clade,Clade>> 
      clades1 = calc1.bipartitions().keySet(),
      clades2 = calc2.bipartitions().keySet();
    double result = 0.0;
    for (UnorderedPair<Clade,Clade> bipart : clades1) if (!clades2.contains(bipart)) result++;
    for (UnorderedPair<Clade,Clade> bipart : clades2) if (!clades1.contains(bipart)) result++;
    return result;
  }
  
  /**
   * All the bipartitions obtained by looping over edges and making a cut.
   * @param tree
   * @return
   */
  public static Set<UnorderedPair<Clade,Clade>> bipartitions(UnrootedTree tree)
  {
    return cladeCalculator(tree, tipIndexer(tree)).bipartitions().keySet();
  }
  
  public static Set<UnorderedPair<Clade,Clade>> bipartitions(UnrootedTree tree, Indexer<TreeNode> tipIndexer)
  {
    return cladeCalculator(tree, tipIndexer).bipartitions().keySet();
  }
  
  public static Indexer<TreeNode> tipIndexer(UnrootedTree p)
  {
    Indexer<TreeNode> tipIndexer = new Indexer<TreeNode>();
    for (TreeNode t : p.leaves())
      tipIndexer.addToIndex(t);
    return tipIndexer;
  }

  private static CladeCalculator cladeCalculator(UnrootedTree ut, Indexer<TreeNode> tipIndexer)
  {
    CladeBuilderFactorGraph fg = new CladeBuilderFactorGraph(ut, tipIndexer);
    SumProduct<TreeNode> sp = new SumProduct<TreeNode>(fg);
    return new CladeCalculator(sp, ut);
  }
  
  /**
   * Cache the bipartitions (used to compute multiple metrics between trees).
   * 
   * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
   */
  private static class CladeCalculator
  {
    private final SumProduct<TreeNode> rawClades;
    private final UnrootedTree ut;
    private Counter<UnorderedPair<Clade,Clade>> bipartitions = null;

    private CladeCalculator(SumProduct<TreeNode> rawClades, UnrootedTree ut)
    {
      this.ut = ut;
      this.rawClades = rawClades;
    }
    
    private Counter<UnorderedPair<Clade,Clade>> bipartitions()
    {
      if (bipartitions != null) return bipartitions;
      bipartitions = new Counter<UnorderedPair<Clade,Clade>>();
      for (UnorderedPair<TreeNode, TreeNode> edge : ut.getTopology().edgeSet())//Graphs.edgeSet(ut.getTopology()))
      {
        double branch = ut.getBranchLength(edge);
        Clade 
          c1 = ((Clade) rawClades.getMessage(edge.getFirst(), edge.getSecond())),
          c2 = ((Clade) rawClades.getMessage(edge.getSecond(), edge.getFirst()));
        if (c1.sortedTips == null && c2.sortedTips == null)
          throw new RuntimeException();
        bipartitions.incrementCount(UnorderedPair.of(c1, c2), branch);
      }
      return bipartitions;
    }
  }
  
  /**
   * A factor graph used internally to efficiently get all the bipartitions clades.
   * 
   * Basically, turn the problem into a sum product algorithm, where pointwise multiplication
   * is union.
   * 
   * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
   *
   */
  private static final class CladeBuilderFactorGraph implements FactorGraph<TreeNode>
  {
    private final UndirectedGraph<TreeNode,?> topology;
    private final Map<TreeNode,UnaryFactor<TreeNode>> map;
    private final Indexer<TreeNode> tipIndexer;
    
    public CladeBuilderFactorGraph(UnrootedTree t, Indexer<TreeNode> indexer)
    {
      this.topology = t.getTopology();
      this.map = Maps.newLinkedHashMap();
      this.tipIndexer = indexer;
      init(t);
    }
    
    private void init(UnrootedTree t)
    {
      for (TreeNode tip : t.leaves())
      {
        int index = tipIndexer.o2i(tip);
        map.put(tip, new Clade(new int[]{index}, tipIndexer));
      }
    }

    @Override
    public UnaryFactor<TreeNode> getUnary(TreeNode node)
    {
      UnaryFactor<TreeNode> result = map.get(node);
      return result;
    }
    @Override public UndirectedGraph<TreeNode,?> getTopology() { return topology; }
    @Override public BinaryFactor<TreeNode> getBinary(TreeNode source, TreeNode destination) { return null; }

    @Override
    public FactorOperations<TreeNode> factorOperations()
    {
      return new FactorOperations<TreeNode>() {

        /**
         * Compute the union. Note: if one of them contains more than half (.sortedTips == null),
         * then taking union will result in a clade with also more than half of the tips (and 
         * therefore where .sortedTips == null.
         * @param otherFactors
         * @return
         */
        @Override
        public UnaryFactor<TreeNode> pointwiseProduct(
            List<? extends UnaryFactor<TreeNode>> otherFactors)
        {
          Indexer<TreeNode> indexer = null;
          int nItems = 0;
          for (UnaryFactor<TreeNode> f : otherFactors)
          {
            Clade other = ((Clade) f);
            if (other.sortedTips == null)
              return other;
            indexer = other.indexer;
            nItems += other.sortedTips.length;
          }
          
          if (nItems > indexer.size()/2)
            return new Clade(null, indexer);
          
          int [] result = new int[nItems];
          int index = 0;
          for (UnaryFactor<TreeNode> f : otherFactors)
            for (int item : ((Clade) f).sortedTips)
            result[index++] = item;

          return new Clade(result, indexer);
        }

        @Override
        public UnaryFactor<TreeNode> marginalize(BinaryFactor<TreeNode> binary,
            List<UnaryFactor<TreeNode>> unariesOnMarginalized)
        {
          return pointwiseProduct(unariesOnMarginalized);
        }
      };
    }
  }
  
  /**
   * The representation of a clade. To save memory, only those with less or equal than
   * half of the tips are maintained. Since we always use these in UnorderedPairs, the one
   * with more than half can be inferred as the complement of the other. Only when
   * the number of tips is even and the bipartition cut it in two halves you will find
   * Unordered pairs where both clades are explicitly represented.
   * 
   * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
   */
  public static class Clade implements UnaryFactor<TreeNode>
  {
    private final int [] sortedTips;
    private final Indexer<TreeNode> indexer;
    public Clade(int[] sortedTips, Indexer<TreeNode> indexer)
    {
      this.indexer = indexer;
      this.sortedTips = sortedTips;
      if (sortedTips != null)
      {
        Arrays.sort(sortedTips);
        for (int i = 0; i < sortedTips.length-1; i++)
          if (sortedTips[i] == sortedTips[i+1])
            throw new RuntimeException();
      }
    }
    @Override
    public String toString()
    {
      if (sortedTips == null)
        return "complement";
      return taxa().toString();
    }
    public Set<TreeNode> taxa()
    {
      Set<TreeNode> result = new HashSet<TreeNode>();
      for (int i : sortedTips)
        result.add(indexer.i2o(i));
      return result;
    }
    @Override
    public int hashCode()
    {
      final int prime = 31;
      int result = 1;
      result = prime * result + Arrays.hashCode(sortedTips);
      return result;
    }
    @Override
    public boolean equals(Object obj)
    {
      if (this == obj)
        return true;
      if (obj == null)
        return false;
      if (getClass() != obj.getClass())
        return false;
      Clade other = (Clade) obj;
      if (other.indexer != this.indexer)
        throw new RuntimeException();
      if (sortedTips == null && other.sortedTips == null)
        return true;
      if (!Arrays.equals(sortedTips, other.sortedTips))
        return false;
      return true;
    }
    @Override
    public double logNormalization()
    {
      return Double.NaN;
    }
  }
}
