package conifer.io;

import conifer.ctmc.expfam.RateMtxNames;
import conifer.ctmc.expfam.SerializedExpFamMixture;
import conifer.models.CNPair;
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


  public static Indexer<String> modelIndexer(final RateMtxNames selectedRateMtx)
  {
    if (selectedRateMtx == null) {
      throw new IllegalArgumentException("model is null!");
    }
    
    return selectedRateMtx.getIndexer();
  }


	public static Indexer<CNPair> CNPairIndexer()
	{
		return CNObservationFactory.defaultFactory().getIndexer();
	}
}


