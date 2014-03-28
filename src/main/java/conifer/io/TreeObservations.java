package conifer.io;

import java.util.List;

import conifer.TreeNode;

/**
 * Maps TreeNodes to data.
 * 
 * Observations can be simple discrete traits, but could also
 * be backed by an alignment to be resampled.
 * 
 * Continuous traits might be another later possibility later on.
 */
public interface TreeObservations
{
  public List<TreeNode> getObservedTreeNodes();

  public Object get(TreeNode leaf);
  public void set(TreeNode leaf, Object data);
  
  public void clear();

  public int nSites();
}