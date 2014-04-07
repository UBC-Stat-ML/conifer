package conifer.ctmc.expfam;

import briefj.Indexer;



public class CTMCStateSpace
{
  /**
   * The partition this space is representing (spaces could
   * be different across different partitions)
   * 
   * Should be null iff there is a single partition in the study.
   */
  public final Object currentPartition;
  
  /**
   * An indexer for the observations (can be partition specific,
   * but is the same for all categories). 
   * 
   * Rationale: partitions could contain different data types 
   * (e.g. points of interaction, points not interacting).
   */
  public final Indexer<String> observationsIndexer;
   
  /**
   * An indexer for the latent states (can be partition specific,
   * but is the same for all categories).
   * 
   * Rationale: partitions could contain different data types 
   * (e.g. points of interaction, points not interacting).
   */
  public final Indexer<?> latentIndexer;
  
  public final int nCategories;
  
  public CTMCStateSpace(
      Indexer<String> observationsIndexer, 
      Indexer<?> latentIndexer,
      int nCategories)
  {
    this(observationsIndexer, latentIndexer, nCategories, null);
  }

  public CTMCStateSpace(
      Indexer<String> observationsIndexer, 
      Indexer<?> latentIndexer,
      int nCategories, 
      Object currentPartition)
  {
    this.currentPartition = currentPartition;
    this.observationsIndexer = observationsIndexer;
    this.latentIndexer = latentIndexer;
    this.nCategories = nCategories;
  }
}