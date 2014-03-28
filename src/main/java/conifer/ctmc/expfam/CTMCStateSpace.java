package conifer.ctmc.expfam;

import briefj.Indexer;



public class CTMCStateSpace
{
  
  public final Object currentPartition;
  
  /**
   * Keys: partitions.
   * 
   * An indexer for the observations (can be partition specific,
   * but is the same for all categories). 
   * 
   * Rationale: partitions could contain different data types 
   * (e.g. points of interaction, points not interacting).
   */
  public final Indexer<String> observationsIndexer;
   
  /**
   * Keys: partitions.
   * 
   * An indexer for the latent states (can be partition specific,
   * but is the same for all categories).
   * 
   * Rationale: partitions could contain different data types 
   * (e.g. points of interaction, points not interacting).
   */
  public final Indexer<Object> latentIndexer;
  
  public final int nCategories;

  public CTMCStateSpace(Object currentPartition,
      Indexer<String> observationsIndexer, Indexer<Object> latentIndexer,
      int nCategories)
  {
    this.currentPartition = currentPartition;
    this.observationsIndexer = observationsIndexer;
    this.latentIndexer = latentIndexer;
    this.nCategories = nCategories;
  }
}