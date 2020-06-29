package conifer

import java.io.File
import java.util.List
import java.util.Map
import java.util.Random
import org.apache.commons.lang3.tuple.Pair
import org.jgrapht.UndirectedGraph
import bayonet.graphs.GraphUtils
import bayonet.marginal.algo.EdgeSorter
import blang.inits.ConstructorArg
import blang.inits.DesignatedConstructor
import blang.mcmc.Samplers
import briefj.collections.Counter
import briefj.collections.UnorderedPair
import com.google.common.collect.Lists
import com.google.common.collect.Maps
import conifer.moves.AllBranchesScaling
import conifer.moves.SPRMove
import conifer.moves.SingleBranchScaling
import conifer.moves.SingleNNI

/** 
 * An unrooted phylogenetic tree.
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 */
@Samplers(SingleNNI, SingleBranchScaling)
class UnrootedTree {
  Map<UnorderedPair<TreeNode, TreeNode>, Double> branchLengths
  UndirectedGraph<TreeNode, UnorderedPair<TreeNode, TreeNode>> topology
  /** 
   * Note if other fields are added, setTo() should be modified as well.
   */

  new(UnrootedTree model) {
    this.topology = GraphUtils::newUndirectedGraph(model.topology)
    this.branchLengths = Maps::newLinkedHashMap(model.branchLengths)
  }

  /** 
   * Create an empty tree.
   */
  new() {
    topology = GraphUtils::newUndirectedGraph()
    branchLengths = Maps::newLinkedHashMap()
  }

  /** 
   * Add a node (species)
   * @param node
   */
  def void addNode(TreeNode node) {
    this.topology.addVertex(node)
  }

  /** 
   * @param node1
   * @param node2
   * @param length
   */
  def void addEdge(TreeNode node1, TreeNode node2, double length) {
    this.topology.addEdge(node1, node2)
    branchLengths.put(new UnorderedPair<TreeNode, TreeNode>(node1, node2), length)
  }

  def void removeEdge(TreeNode n1, TreeNode n2) {
    topology.removeEdge(n1, n2)
    branchLengths.remove(UnorderedPair::of(n1, n2))
  }

  def void removeEdge(UnorderedPair<TreeNode, TreeNode> edge) {
    removeEdge(edge.getFirst(), edge.getSecond())
  }

  /** 
   * @return
   */
  def UndirectedGraph<TreeNode, UnorderedPair<TreeNode, TreeNode>> getTopology() {
    return topology
  }

  /** 
   * @return
   */
  def Map<UnorderedPair<TreeNode, TreeNode>, Double> getBranchLengths() {
    return branchLengths
  }

  /** 
   * Return the branch corresponding to the unordered pair {node1, node2}
   * @param node1
   * @param node2
   * @return
   */
  def Double getBranchLength(TreeNode node1, TreeNode node2) {
    return branchLengths.get(new UnorderedPair<TreeNode, TreeNode>(node1, node2))
  }

  def double getBranchLength(UnorderedPair<TreeNode, TreeNode> edge) {
    return branchLengths.get(edge)
  }

  /** 
   * Modifies in place the length of a branch.
   * @param edge
   * @param newValue
   */
  def void updateBranchLength(UnorderedPair<TreeNode, TreeNode> edge, double newValue) {
    if(!branchLengths.containsKey(edge)) throw new RuntimeException();
    branchLengths.put(edge, newValue)
  }

  /** 
   * Orient the edges relative to the provided root.
   * For each pair p, p.getLeft() contains the parent,
   * and p.getLeft() contains the children.
   * @param root
   * @return The oriented edges
   */
  def List<Pair<TreeNode, TreeNode>> getRootedEdges(TreeNode root) {
    return EdgeSorter::newEdgeSorter(getTopology(), root).backwardMessages()
  }

  /** 
   * Performs a Nearest neighbor interchange move.
   * @param moved1
   * @param fixed1
   * @param moved2
   * @param fixed2
   */
  def void interchange(TreeNode moved1, TreeNode fixed1, TreeNode moved2, TreeNode fixed2) {
    modifyBranch(fixed1, moved1, fixed2)
    modifyBranch(fixed2, moved2, fixed1)
  }

  def private void modifyBranch(TreeNode oldFixed, TreeNode moved, TreeNode newFixed) {
    var double branchLength = getBranchLength(oldFixed, moved)
    removeEdge(oldFixed, moved)
    addEdge(moved, newFixed, branchLength)
  }

  /** 
   * Set the value of the topology and branch length of this tree to those of the provided
   * tree (no copy performed, just two references set)
   * @param otherTree
   */
  def void setTo(UnrootedTree otherTree) {
    this.topology = otherTree.topology
    this.branchLengths = otherTree.branchLengths
  }

  /** 
   * See UnrootedTreeUtils.leaves()
   * @return
   */
  def List<TreeNode> leaves() {
    return GraphUtils::leaves(getTopology())
  }

  /** 
   * See UnrootedTreeUtils.allTotalBranchLengthDistances()
   * @return
   */
  def Counter<UnorderedPair<TreeNode, TreeNode>> allTotalBranchLengthDistances() {
    return UnrootedTreeUtils::allTotalBranchLengthDistances(this)
  }

  /** 
   * Reads the contents of the given file and parse it as a newick string.
   * See {@link #fromNewickString(String) fromNewickStrin(String)}
   */
  @DesignatedConstructor def static UnrootedTree fromNewick(@ConstructorArg("newickFile") File f) {
    return UnrootedTreeUtils::fromNewick(f)
  }

  /** 
   * Parse a tree from a string containing a newick specification.
   * Limitations:
   * - Names of taxa should have no space, start with ["a"-"z","A"-"Z","_"], and have character from ["a"-"z","A"-"Z","_","-","0"-"9","."] for the rest
   * - Leaves should be named (not required for internal)
   * - Leaf names should be distrinct
   * - All edges should have a branch length attached to it
   * For example: UnrootedTree.fromNewickString("((A:1.0,Z:2.0):3.0,(B:4.0,C:5.0):6.0,X:100);");
   * @param string
   * @return
   */
  def static UnrootedTree fromNewickString(String string) {
    return UnrootedTreeUtils::fromNewickString(string)
  }

  def String toNewick() {
    return UnrootedTreeUtils::toNewick(this)
  }

  override String toString() {
    return toNewick()
  }

  /** 
   * Iterate the edge (oriented with the provided root) and add a dummy internal node on 
   * each edge, except for edges connected to current. 
   * The nodes are placed at a uniform fractions from the bottom node of each edge. 
   * This modifies the tree in place.
   */
  def List<TreeNode> addAuxiliaryInternalNodes(Random rand, TreeNode current) {
    var List<TreeNode> result = Lists::newArrayList()
    for (Pair<TreeNode, TreeNode> edge : getRootedEdges(current)) {
      if(edge.getLeft().equals(current) || edge.getRight().equals(current))
        result.add(current)
      else {
        var double ratio = rand.nextDouble()
        var double originalBL = getBranchLength(edge.getLeft(), edge.getRight())
        var double bottomBL = ratio * originalBL
        var double top_BL = (1.0 - ratio) * originalBL
        removeEdge(edge.getLeft(), edge.getRight())
        var TreeNode dummyNode = TreeNode::nextUnlabelled()
        addNode(dummyNode)
        // left = complete top
        addEdge(edge.getLeft(), dummyNode, top_BL)
        addEdge(dummyNode, edge.getRight(), bottomBL)
        result.add(dummyNode)
      }
    }
    return result
  }

  /** 
   * This assumes that removedRoot has exactly 2 neighbors, n1 and n2 (hence they form a chain ..n1--removedRoot--n2..
   * Simplify the tree in place into ..n1--n2.. with the new branch length equal to the
   * sum of the two removed edges.
   * @param n1
   * @param removedRoot
   * @param n2
   */
  def void simplify(TreeNode n1, TreeNode removedRoot, TreeNode n2) {
    if(topology.edgesOf(removedRoot).size() !== 2) throw new RuntimeException(
      "Simplify assumes that the node to be removed has exactly 2 neighbors.");
    var double branchSum = getBranchLength(n1, removedRoot) + getBranchLength(removedRoot, n2)
    removeEdge(n1, removedRoot)
    removeEdge(removedRoot, n2)
    topology.removeVertex(removedRoot)
    addEdge(n1, n2, branchSum)
  }

  /** 
   * This splits the tree into two parts relative to the edge e=(removedRoot, detached).
   * This will return a subtree containing e and the subtree on the side of the node detached
   * relative to e. 
   * This instance will be modified in place to remove all the edges and nodes in the returned
   * tree, except for the node removedRoot.
   */
  def UnrootedTree prune(TreeNode removedRoot, TreeNode detachedNode) {
    var double branchLen = getBranchLength(removedRoot, detachedNode)
    removeEdge(removedRoot, detachedNode)
    var UnrootedTree result = new UnrootedTree()
    result.addNode(detachedNode)
    for (Pair<TreeNode, TreeNode> orientedEdge : getRootedEdges(detachedNode)) {
      result.addNode(orientedEdge.getRight())
      result.addEdge(orientedEdge.getLeft(), orientedEdge.getRight(),
        getBranchLength(orientedEdge.getLeft(), orientedEdge.getRight()))
      this.removeEdge(orientedEdge.getLeft(), orientedEdge.getRight())
      this.topology.removeVertex(orientedEdge.getRight())
    }
    this.topology.removeVertex(detachedNode)
    result.addNode(removedRoot)
    result.addEdge(detachedNode, removedRoot, branchLen)
    return result
  }

  /** 
   * Take the subtree prunedSubtree rooted at prunedSubtreeRoot,
   * and add it in place, essentially merging node specified by variable attachment with node
   * specified by prunedSubtreeRoot (keeping the name in variable attachment)
   * @param prunedSubtree
   * @param prunedSubtreeRoot
   * @param attachment
   */
  def void regraft(UnrootedTree prunedSubtree, TreeNode prunedSubtreeRoot, TreeNode attachment) {
    for (Pair<TreeNode, TreeNode> orientedEdge : prunedSubtree.getRootedEdges(prunedSubtreeRoot)) {
      var TreeNode topNode = orientedEdge.getLeft()
      var TreeNode botNode = orientedEdge.getRight()
      var double branchLength = prunedSubtree.getBranchLength(topNode, botNode)
      if(topNode.equals(prunedSubtreeRoot)) topNode = attachment
      this.addNode(botNode)
      this.addEdge(topNode, botNode, branchLength)
    }
  }

  /** 
   * Remove all internal nodes with exactly two neighbors, changing edges and 
   * branch lengths accordingly to keep the same interpretation of the tree
   * (under reversible models)
   */
  def void simplify() {
    UnrootedTreeUtils::simplify(this)
  }
}
