package conifer.factors;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.math.SamplingUtils;
import blang.core.RealDistribution;
import briefj.collections.UnorderedPair;

import com.google.common.collect.Lists;

import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;



public class NonClockTreePriorUtils<P>
{
  /**
   * Generate a tree with a topology uniformly distributed
   * over bifurcating unrooted tree, using:
   * 
   *   The generation of random, binary unordered trees
   *   George W. Furnas
   *     see 2.1.2 (p.204)
   */
  public static UnrootedTree sample(
      Random random, 
      RealDistribution branchDistribution,
      Collection<TreeNode> leaves)
  {
    UnrootedTree result = new UnrootedTree();
    
    List<TreeNode> shuffled = Lists.newArrayList(leaves);
    Collections.shuffle(shuffled, random);
    Queue<TreeNode> queue = Lists.newLinkedList(shuffled);
    
    if (queue.isEmpty())
      return result;
    
    TreeNode leaf1 = queue.poll();
    result.addNode(leaf1);
    
    if (queue.isEmpty())
      return result;
    
    TreeNode leaf2 = queue.poll();
    result.addNode(leaf2);
    result.addEdge(leaf1, leaf2, Double.NaN);
    
    while (!queue.isEmpty())
    {
      // pick a random edge
      UnorderedPair<TreeNode, TreeNode> edge = SamplingUtils.uniformFromCollection(random, result.getTopology().edgeSet());
      TreeNode internal = TreeNode.nextUnlabelled();
      TreeNode newLeaf = queue.poll();
      result.removeEdge(edge);
      result.addNode(newLeaf);
      result.addNode(internal);
      result.addEdge(newLeaf, internal, Double.NaN);
      result.addEdge(internal, edge.getFirst(), Double.NaN);
      result.addEdge(internal, edge.getSecond(), Double.NaN);
    }
    
    for (UnorderedPair<TreeNode, TreeNode> edge : result.getTopology().edgeSet())
      result.updateBranchLength(edge, branchDistribution.sample(random));
    
    return result;
  }
  
  public static UnrootedTree sampleBalancedUnrootedBinaryTree(
	      double defaultLength,
	      Collection<TreeNode> leaves)
	  {
	    UnrootedTree result = new UnrootedTree();
	    
	    List<TreeNode> shuffled = Lists.newArrayList(leaves);
	    Queue<TreeNode> queue = Lists.newLinkedList(shuffled);
	   
	    if (queue.isEmpty())
	      return result;
	    
	    int nLeaves = leaves.size();
	    double initialHeight = Math.floor(Math.log(nLeaves)/Math.log(2));
	    int nBinaryLeaves = (int) Math.pow(2, initialHeight);
	    int  remainder = nLeaves - nBinaryLeaves;
	    
	    List<TreeNode> firstNRemainderLeaves = shuffled.subList(0, (remainder));
	    List<TreeNode> binaryLeaves = shuffled.subList(remainder,(shuffled.size()));
	    Queue<TreeNode> queueRemainder = Lists.newLinkedList(firstNRemainderLeaves);
	    Queue<TreeNode> binaryQueue = Lists.newLinkedList(binaryLeaves);
	   
	    // use Recursion algorithm to build a complete binary tree	  
	    Queue<TreeNode> parentQueue = binaryQueue;
	   
	    // make a copy of a Queue parentQueue
	    Queue<TreeNode> parentQueueCopy = Lists.newLinkedList(binaryLeaves);
	    
	    Queue<TreeNode> childQueue = Lists.newLinkedList();
	    int count = 0;
	    Map<TreeNode, TreeNode> leavesEdgeSet = new HashMap<TreeNode, TreeNode>(); 
	    
	    while(parentQueue.size()>1){   
	    	TreeNode leaf1 = parentQueue.poll();
	    	result.addNode(leaf1);
	 
	    	TreeNode leaf2 = parentQueue.poll();
	    	result.addNode(leaf2);
	    	TreeNode internal = TreeNode.nextUnlabelled();
	    	childQueue.add(internal);
	    	result.addNode(internal);
	    	
	    	result.addEdge(internal, leaf1, defaultLength);
	    	result.addEdge(internal, leaf2, defaultLength);
	    	
	    	if(count < nBinaryLeaves/2 ){
	    		leavesEdgeSet.put(leaf1, internal);
	    		leavesEdgeSet.put(leaf2, internal);
	    		count++;    		
	    	}
	    	
	    	if(parentQueue.isEmpty()){
	    		
	    		parentQueue = Lists.newLinkedList(childQueue);
	    		
	    		if(childQueue.size() == nBinaryLeaves/2){
	    			
	    			// for each leaf in the remainderQueue, we randomly choose one of the leaves in the current binary tree
	    			// break the current branch with the chosen leaf and its parent-the old internal node, 
	    			// then we add one new internal node, then we build three new branches, the branch connecting the new 
	    			// internal node and the leaf in the remainderQueue, a new branch connecting the new internal node and the 
	    			// chosen leaf, and a branch connecting the new internal node and the old internal node
	    			while(!queueRemainder.isEmpty()){
	    				
	    		    		TreeNode leaf3 = queueRemainder.poll();
	    		    		result.addNode(leaf3);
	    		    	
	    		    		// create a new internal node
	    		    		TreeNode internalNew = TreeNode.nextUnlabelled();
	    		    		result.addNode(internalNew);
	    		    		// build a new branch between the new internal node and the leaf in the queueRemainder
	    		    		result.addEdge(leaf3, internalNew, defaultLength);
	    		    	    
	    		    		// shuffle the leaves in the binary tree randomly
	    		    		// randomly pick a leaf from the binary tree
	    		    		// once a leaf has been chosen, poll it out from the Queue
	    		    		TreeNode binaryLeaf = parentQueueCopy.poll();
	    		    		// build an edge between the new internal node and the binaryLeaf
	    		    		result.addEdge(binaryLeaf, internalNew, defaultLength);
	    		    		 		    		
	    		    		// get the parent internal node of the chosen binaryLeaf
	    		    		// build an edge between the parent internal node and the new internal node
	    		    		
	    		    	    TreeNode internalOld = leavesEdgeSet.get(binaryLeaf);
	    		    	    UnorderedPair<TreeNode, TreeNode> edge = new UnorderedPair<TreeNode, TreeNode>(binaryLeaf, internalOld);
	    		    	    result.removeEdge(edge);
	    		    	    result.addNode(internalOld);
	    		    	    result.addEdge(internalOld, internalNew, defaultLength);
	    		    	    
	    		    	    // if we insert a new internal node for this leaf, the original leaf reference should be updated to the new internal node
	    		    	    if(leaf1.equals(binaryLeaf)){
	    		    	    		
	    		    	    		leaf1 = internalNew;
	    		    	    		
	    		    	    }
	    		    	    
	    		    	    if(leaf2.equals(binaryLeaf)){
    		    	    		
	    		    	    		leaf2 = internalNew;
    		    	    		
	    		    	    }
	    		    	    
	    		    }
	    			
	    		}
	    		if(childQueue.size()>1){
	    			childQueue.clear();
	    		}else{
	    			// remove the root of the binary tree
	    			//result.simplify(leaf1, internal, leaf2);
	    			result.removeEdge(internal, leaf1);
	    			result.removeEdge(internal, leaf2);
	    			result.getTopology().removeVertex(internal);
	    			result.addEdge(leaf1, leaf2, defaultLength);
	       		}
	    		
	    	}
	    	    		    	
	    }
	    return result;
	 
	  }
  
  
  public static <T> Pair<T,T> popRandomPair(Random rand, List<T> items)
  {
    if (items.size() < 2) 
      throw new RuntimeException();
    List<Integer> indices = sampleWithoutReplacement(rand, items.size(), 2);
    Collections.sort(indices);
    Pair<T,T> result = Pair.of(items.get(indices.get(0)), items.get(indices.get(1)));
    items.remove((int)indices.get(1)); // remove second before first to avoid shifts
    items.remove((int)indices.get(0));
    return result;
  }
  
  /**
   * Returns n samples from {0, ..., s-1}, without replacement
   * 
   * TODO: could be more efficient
   * TODO: move to Bayonet
   * 
   * @param n
   * @param s
   * @return
   */
  public static List<Integer> sampleWithoutReplacement(Random rand, int s, int n)
  {
    if (n > s || s < 0 || n < 0)
      throw new RuntimeException();
    List<Integer> list = new ArrayList<Integer>(s),
                  result=new ArrayList<Integer>(n);
    for (int i = 0; i < s; i++)
      list.add(i);
    Collections.shuffle(list,rand);
    for (int i = 0; i < n; i++)
      result.add(list.get(i));
    return result;
  }

public P branchDistributionParameters;

public static void main(String[] args){
	
	// create a tree of 18 leaves and see if the generated tree is correct or not
	// create a collection of leaves
	Random rand = new Random(1);
	double defaultLength = 1.0;
	Collection<TreeNode> leaves = TopologyUtils.syntheticTaxaList((int) 15);
	UnrootedTree unrootedTree = sampleBalancedUnrootedBinaryTree(rand, defaultLength,leaves);
	System.out.println(unrootedTree.toNewick());	
	}

}
