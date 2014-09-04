package conifer.io;

import java.util.ArrayList;
import java.util.List;

import conifer.models.CNPair;
import briefj.Indexer;

public class CNObservationFactory {

	private transient Indexer<CNPair> _indexer;
	private final List<CNPair> orderedSymbols = allPossibleCNPairs();
	
	public static final int maximumNormalCopyNumber = 6;
	public static final int maximumMutantCopyNumber = 6;

	public Indexer<CNPair> getIndexer()
	{
		if (_indexer == null) 
			_indexer = new Indexer<CNPair>(orderedSymbols);
		return _indexer;
	}

	public static List<CNPair> allPossibleCNPairs()
	{	
		List<CNPair> result = new ArrayList<CNPair>();

		for (int i = 0; i < maximumNormalCopyNumber; i++) {
			for (int j = 0; j < maximumMutantCopyNumber; j++) {
				result.add(CNPair.withCounts(i, j));
			}
		}

		return result;
	}
	
	public static CNObservationFactory defaultFactory() {
		return new CNObservationFactory();
	}

  public double[][] site2CharacterIndicators(String rawString)
  {
    // TODO: needed for loading observations 
    return null;
  }
}
