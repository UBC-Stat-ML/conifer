package conifer.io;

import conifer.models.CNPair;
import briefj.Indexer;

public class Indexers {
	public static Indexer<String> dnaIndexer() {
		return PhylogeneticObservationFactory.nucleotidesFactory().getIndexer();
	}

	public static Indexer<String> proteinIndexer() {
		return PhylogeneticObservationFactory.proteinFactory().getIndexer();
	}

	public static Indexer<String> proteinPairIndexer() {
		return PhylogeneticObservationFactory.proteinPairFactory().getIndexer();
	}

	public static Indexer<String> copyNumberCTMCIndexer() {
		return PhylogeneticObservationFactory.copyNumberCTMCFactory().getIndexer();
	}

	public static Indexer<String> copyNumberEmissionIndexer() {
		return PhylogeneticObservationFactory.copyNumberEmissionFactory().getIndexer();
	}

	public static Indexer<CNPair> CNPairIndexer() {
		return CNObservationFactory.defaultFactory().getIndexer();
	}

	public static void main(String[] args) {
		Indexer<String> i = Indexers.copyNumberCTMCIndexer();
		for (Object j : i.objectsList()) {
			System.out.println(j.toString());
		}
	}

}
