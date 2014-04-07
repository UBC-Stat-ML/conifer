package conifer.io;

import briefj.Indexer;



public class Indexers
{
  public static Indexer<String> dnaIndexer()
  {
    return PhylogeneticObservationFactory.nucleotidesFactory().getIndexer();
  }
}
