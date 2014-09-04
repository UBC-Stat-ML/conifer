package conifer.io;

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

	public static Indexer<String> proteinPairIndexer()
	{
		return PhylogeneticObservationFactory.proteinPairFactory().getIndexer();
	}

	public static Indexer<CNPair> CNPairIndexer()
	{
		return CNObservationFactory.defaultFactory().getIndexer();
	}
}


