package conifer.ctmc.expfam;



public class CTMCState
{
  
  /**
   * A category is a connected component in the CTMC.
   * 
   * E.g.: site-specific rate variation.
   */
  public final int categoryIndex;
  
  /**
   * Any other latent information that do not follow
   * connected component-ness.
   * 
   * E.g.: there would be eight of these in a covarion,
   * four of these in a GTR+Gamma2
   * model.
   */
  public final Object latent;
  
  /**
   * A partition is a subset of the dataset.
   * 
   * Should be null iff there is a single partition.
   * 
   * E.g.: Which gene this state belongs to. 
   * Or position index in a codon (e.g. for a CAT model).
   */
  public final Object partition;
  
  
  public CTMCState(
      int categoryIndex, 
      Object latent,
      Object partition)
  {
    this.categoryIndex = categoryIndex;
    this.latent = latent;
    this.partition = partition;
  }
  
  @Override
  public String toString()
  {
    return latent + "(" + categoryIndex + (partition == null ? "" : "," + partition) + ")";
  }

  @Override
  public int hashCode()
  {
    final int prime = 31;
    int result = 1;
    result = prime * result + categoryIndex;
    result = prime * result + ((latent == null) ? 0 : latent.hashCode());
    result = prime * result + ((partition == null) ? 0 : partition.hashCode());
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
    CTMCState other = (CTMCState) obj;
    if (categoryIndex != other.categoryIndex)
      return false;
    if (latent == null)
    {
      if (other.latent != null)
        return false;
    } else if (!latent.equals(other.latent))
      return false;
    if (partition == null)
    {
      if (other.partition != null)
        return false;
    } else if (!partition.equals(other.partition))
      return false;
    return true;
  }
  
  
 
}