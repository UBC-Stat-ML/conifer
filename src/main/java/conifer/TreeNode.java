package conifer;


/**
 * The name of a node on a tree. Can be an informative string for the leaves of the 
 * tree, or some id for the internal nodes.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public final class TreeNode
{
  public static TreeNode withLabel(String name)
  {
    return new TreeNode(name);
  }
  
  /**
   * Creates a unique internal node.
   * @return
   */
  public static TreeNode nextUnlabelled()
  {
    synchronized (INTERNAL_PREFIX)
    {
      return withLabel(INTERNAL_PREFIX + (nextId++));
    }
  }
  
  @Override
  public String toString()
  {
    return description;
  }
  
  public boolean isLabelled() { return isLabelled; }
  
  private static int nextId = 0;
  
  // Note: no need for toString, hashCode because strings are interned
  
  private final String description;
  private final boolean isLabelled;
  
  private static final String INTERNAL_PREFIX = "unlabelled_";
  
  private TreeNode(String description)
  {
    this.description = description.intern();
    this.isLabelled = !description.startsWith(INTERNAL_PREFIX);
  }
}
