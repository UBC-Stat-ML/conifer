package conifer.io;

import conifer.ctmc.expfam.RateMtxNames;
import conifer.ctmc.expfam.SerializedExpFamMixture;
import briefj.Indexer;



public class Indexers
{
  public static Indexer<String> dnaIndexer()
  {
    return PhylogeneticObservationFactory.nucleotidesFactory().getIndexer();
  }
  public static Indexer<String> proteinIndexer()
  {
    return PhylogeneticObservationFactory.proteinFactory().getIndexer();
  }

  public static Indexer<String> proteinPairIndexer()
  {
    return PhylogeneticObservationFactory.proteinPairFactory().getIndexer();
  }

  public static Indexer<String> modelIndexer(final String model)
  {
    if (model == null) {
      throw new IllegalArgumentException("model is null!");
    }

    final RateMtxNames something = RateMtxNames.fromString(model);

    if (something == null) {
      return new Indexer();
    }

    switch(something) {
    case KIMURA1980:
      return Indexers.dnaIndexer();
    case ACCORDANCE:
      return Indexers.proteinIndexer();
    case PAIR:
      return Indexers.proteinPairIndexer();
    case POLARITY:
      return Indexers.proteinIndexer();
    case POLARITYSIZE:
      return Indexers.proteinIndexer();

    default:
      return new Indexer();    
    }   

  }

}


