package conifer;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;
import org.jgrapht.UndirectedGraph;

import bayonet.graphs.GraphUtils;
import bayonet.marginal.algo.EdgeSorter;
import blang.annotations.Processors;
import blang.annotations.Samplers;
import briefj.BriefCollections;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import conifer.moves.AllBranchesScaling;
import conifer.moves.SPRMove;
import conifer.moves.SingleBranchScaling;
import conifer.moves.SingleNNI;
import conifer.processors.TotalTreeLengthProcessor;
import conifer.processors.TreeDiameterProcessor;
import conifer.processors.TreeDistanceProcessor;


/**
 * An unrooted phylogenetic tree.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
@Samplers({
  SingleNNI.class, 				// only topology
  SingleBranchScaling.class, // only branch length
  AllBranchesScaling.class, // affect only branch length
  SPRMove.class				// only topology
})
@Processors({
  TotalTreeLengthProcessor.class,
  TreeDiameterProcessor.class,
  TreeDistanceProcessor.class
})
public class UnrootedTree
{
  /**
   * Note if other fields are added, setTo() should be modified as well.
   */
  private UndirectedGraph<TreeNode, UnorderedPair<TreeNode, TreeNode>> topology;
  private Map<UnorderedPair<TreeNode, TreeNode>,Double> branchLengths;
  
  public UnrootedTree(UnrootedTree model)
  {
    this.topology = GraphUtils.newUndirectedGraph(model.topology);
    this.branchLengths = Maps.newLinkedHashMap(model.branchLengths);
  }
  
  /**
   * Create an empty tree.
   */
  public UnrootedTree()
  {
    topology = GraphUtils.newUndirectedGraph();
    branchLengths = Maps.newLinkedHashMap();
  }
  
  /**
   * Add a node (species)
   * 
   * @param node
   */
  public void addNode(TreeNode node)
  {
    this.topology.addVertex(node);
  }
  
  /**
   * @param node1
   * @param node2
   * @param length
   */
  public void addEdge(TreeNode node1, TreeNode node2, double length)
  {
    this.topology.addEdge(node1, node2);
    branchLengths.put(new UnorderedPair<TreeNode, TreeNode>(node1, node2), length);
  }
  
  public void removeEdge(TreeNode n1, TreeNode n2)
  {
    topology.removeEdge(n1, n2);
    branchLengths.remove(UnorderedPair.of(n1, n2));
  }
  
  public void removeEdge(UnorderedPair<TreeNode, TreeNode> edge)
  {
    removeEdge(edge.getFirst(), edge.getSecond());
  }
  
  /**
   * 
   * @return
   */
  public UndirectedGraph<TreeNode, UnorderedPair<TreeNode, TreeNode>> getTopology()
  {
    return topology;
  }
  
  /**
   * 
   * @return
   */
  public Map<UnorderedPair<TreeNode, TreeNode>, Double> getBranchLengths()
  {
    return branchLengths;
  }
  
  /**
   * Return the branch corresponding to the unordered pair {node1, node2}
   * @param node1
   * @param node2
   * @return
   */
  public Double getBranchLength(TreeNode node1, TreeNode node2)
  {
    return branchLengths.get(new UnorderedPair<TreeNode, TreeNode>(node1, node2));
  }
  
  public double getBranchLength(UnorderedPair<TreeNode, TreeNode> edge)
  {
    return branchLengths.get(edge);
  }
 
  
  /**
   * Modifies in place the length of a branch.
   * @param edge
   * @param newValue
   */
  public void updateBranchLength(UnorderedPair<TreeNode, TreeNode> edge,
      double newValue)
  {
    if (!branchLengths.containsKey(edge))
      throw new RuntimeException();
    branchLengths.put(edge, newValue);
  }
  
  /**
   * Orient the edges relative to the provided root.
   * 
   * For each pair p, p.getLeft() contains the parent,
   * and p.getLeft() contains the children.
   * 
   * @param root
   * @return The oriented edges
   */
  public List<Pair<TreeNode,TreeNode>> getRootedEdges(TreeNode root)
  {
    return EdgeSorter.newEdgeSorter(getTopology(), root).backwardMessages();
  }
  
  /**
   * Performs a Nearest neighbor interchange move.
   * 
   * @param moved1
   * @param fixed1
   * @param moved2
   * @param fixed2
   */
  public void interchange(TreeNode moved1, TreeNode fixed1, TreeNode moved2, TreeNode fixed2)
  {
    modifyBranch(fixed1, moved1, fixed2);
    modifyBranch(fixed2, moved2, fixed1);
  }
  
  public TreeNode getTreeNode(String label)
  {
      if (topology.containsVertex(TreeNode.withLabel(label)))
          for(TreeNode n : topology.vertexSet())
          {
              if (n.toString() == "root")
                  return n;
          }
      return null;
  }
  
  public TreeNode getInternalNode()
  {
      List<TreeNode> internalNodes = GraphUtils.internalNodes(topology);
      return BriefCollections.pick(internalNodes);
  }
  
  private void modifyBranch(
      TreeNode oldFixed, 
      TreeNode moved, 
      TreeNode newFixed)
  {
    double branchLength = getBranchLength(oldFixed, moved);
    removeEdge(oldFixed, moved);
    addEdge(moved, newFixed, branchLength);
  }
  
  /**
   * Set the value of the topology and branch length of this tree to those of the provided
   * tree (no copy performed, just two references set)
   * @param otherTree
   */
  public void setTo(UnrootedTree otherTree)
  {
    this.topology = otherTree.topology;
    this.branchLengths = otherTree.branchLengths;
  }

  /**
   * See UnrootedTreeUtils.leaves()
   * @return
   */
  public List<TreeNode> leaves()
  {
    return GraphUtils.leaves(getTopology());
  }
  
  /**
   * See UnrootedTreeUtils.allTotalBranchLengthDistances()
   * @return
   */
  public Counter<UnorderedPair<TreeNode,TreeNode>> allTotalBranchLengthDistances()
  {
    return UnrootedTreeUtils.allTotalBranchLengthDistances(this);
  }
  
  /**
   * Reads the contents of the given file and parse it as a newick string.
   * 
   * See {@link #fromNewickString(String) fromNewickStrin(String)}
   */
  public static UnrootedTree fromNewick(File f)
  {
    return UnrootedTreeUtils.fromNewick(f);
  }
  
  /**
   * Parse a tree from a string containing a newick specification.
   * 
   * Limitations:
   * - Names of taxa should have no space, start with ["a"-"z","A"-"Z","_"], and have character from ["a"-"z","A"-"Z","_","-","0"-"9","."] for the rest
   * - Leaves should be named (not required for internal)
   * - Leaf names should be distrinct
   * - All edges should have a branch length attached to it
   * 
   * For example: UnrootedTree.fromNewickString("((A:1.0,Z:2.0):3.0,(B:4.0,C:5.0):6.0,X:100);");
   * @param string
   * @return
   */
  public static UnrootedTree fromNewickString(String string)
  {
    return UnrootedTreeUtils.fromNewickString(string);
  }
  
  public String toNewick()
  {
    return UnrootedTreeUtils.toNewick(this);
  }
  
  @Override
  public String toString()
  {
    return toNewick();
  }

  /**
   * Iterate the edge (oriented with the provided root) and add a dummy internal node on 
   * each edge, except for edges connected to current. 
   * 
   * The nodes are placed at a uniform fractions from the bottom node of each edge. 
   * This modifies the tree in place.
   */
  public List<TreeNode> addAuxiliaryInternalNodes(Random rand, TreeNode current)
  {
    List<TreeNode> result = Lists.newArrayList();
    
    for (Pair<TreeNode,TreeNode> edge : getRootedEdges(current))
    {
      if (edge.getLeft().equals(current) || edge.getRight().equals(current))
        result.add(current);
      else
      {
        double ratio = rand.nextDouble();
        double originalBL = getBranchLength(edge.getLeft(), edge.getRight());
        double
          bottomBL = ratio * originalBL,
          top_BL = (1.0 - ratio) * originalBL;
  
        removeEdge(edge.getLeft(), edge.getRight());
        TreeNode 
          dummyNode = TreeNode.nextUnlabelled();
        addNode(dummyNode);
        // left = complete top
        addEdge(edge.getLeft(), dummyNode, top_BL);
        addEdge(dummyNode, edge.getRight(), bottomBL);
        result.add(dummyNode);
      }
    }
    
    return result;
  }

  /**
   * This assumes that removedRoot has exactly 2 neighbors, n1 and n2 (hence they form a chain ..n1--removedRoot--n2..
   * 
   * Simplify the tree in place into ..n1--n2.. with the new branch length equal to the
   * sum of the two removed edges.
   * 
   * @param n1
   * @param removedRoot
   * @param n2
   */
  public void simplify(TreeNode n1, TreeNode removedRoot, TreeNode n2)
  {
    if (topology.edgesOf(removedRoot).size() != 2)
      throw new RuntimeException("Simplify assumes that the node to be removed has exactly 2 neighbors.");
    double branchSum = getBranchLength(n1, removedRoot) + getBranchLength(removedRoot, n2);
    removeEdge(n1, removedRoot);
    removeEdge(removedRoot, n2);
    topology.removeVertex(removedRoot);
    addEdge(n1, n2, branchSum);
  }

  /**
   * This splits the tree into two parts relative to the edge e=(removedRoot, detached).
   * 
   * This will return a subtree containing e and the subtree on the side of the node detached
   * relative to e. 
   * 
   * This instance will be modified in place to remove all the edges and nodes in the returned
   * tree, except for the node removedRoot.
   * 
   */
  public UnrootedTree prune(TreeNode removedRoot, TreeNode detachedNode)
  {
    double branchLen = getBranchLength(removedRoot, detachedNode);
    removeEdge(removedRoot, detachedNode);
    
    UnrootedTree result = new UnrootedTree();

    result.addNode(detachedNode);
    
    for (Pair<TreeNode,TreeNode> orientedEdge : getRootedEdges(detachedNode))
    {
      result.addNode(orientedEdge.getRight());
      result.addEdge(orientedEdge.getLeft(), orientedEdge.getRight(), getBranchLength(orientedEdge.getLeft(), orientedEdge.getRight()));
      this.removeEdge(orientedEdge.getLeft(), orientedEdge.getRight());
      this.topology.removeVertex(orientedEdge.getRight());
    }
    
    this.topology.removeVertex(detachedNode);
    result.addNode(removedRoot);
    result.addEdge(detachedNode, removedRoot, branchLen);
    
    return result;
  }

  /**
   * Take the subtree prunedSubtree rooted at prunedSubtreeRoot,
   * and add it in place, essentially merging node specified by variable attachment with node
   * specified by prunedSubtreeRoot (keeping the name in variable attachment)
   * @param prunedSubtree
   * @param prunedSubtreeRoot
   * @param attachment
   */
  public void regraft(UnrootedTree prunedSubtree, TreeNode prunedSubtreeRoot,
      TreeNode attachment)
  {
    for (Pair<TreeNode,TreeNode> orientedEdge : prunedSubtree.getRootedEdges(prunedSubtreeRoot))
    {
      TreeNode 
        topNode = orientedEdge.getLeft(),
        botNode = orientedEdge.getRight();
      double branchLength = prunedSubtree.getBranchLength(topNode, botNode);
      if (topNode.equals(prunedSubtreeRoot))
        topNode = attachment;
      this.addNode(botNode);
      this.addEdge(topNode, botNode, branchLength);
    }
  }

  /**
   * Remove all internal nodes with exactly two neighbors, changing edges and 
   * branch lengths accordingly to keep the same interpretation of the tree
   * (under reversible models)
   */
  public void simplify()
  {
    UnrootedTreeUtils.simplify(this);
  }




}
