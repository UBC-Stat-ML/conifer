package conifer.io;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.ejml.simple.SimpleMatrix;


import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.gson.Gson;

import bayonet.math.EJMLUtils;
import briefj.BriefCollections;
import briefj.BriefIO;
import briefj.BriefLog;
import briefj.Indexer;



public class PhylogeneticObservationFactory
{
  private final List<String> orderedSymbols;
  private final Map<String, Set<String>> ambiguousSymbols;
  private final boolean caseSensitive;
  
  private transient Integer chunkLength = null;
  private transient Indexer<String> _indexer;
  private transient Map<String,double[]> _indicators;
  private transient double[] _unknownIndicator;
  
  public double [][] site2CharacterIndicators(String sequence)
  {
    Map<String,double[]> indicators = getIndicators();
    if (sequence.length() % chunkLength != 0)
      throw new RuntimeException("Sequence length was expected to be a multiple of " + chunkLength);
    int nChunks = sequence.length() / chunkLength;
    double [][] result = new double[nChunks][];
    for (int chunkIndex = 0; chunkIndex < nChunks; chunkIndex++)
    {
      final int startPoint = chunkIndex * chunkLength;
      String currentChunk = sequence.substring(startPoint, startPoint + chunkLength);
      if (!caseSensitive)
        currentChunk = currentChunk.toUpperCase();
      double [] currentIndicator = indicators.get(currentChunk);
      if (currentIndicator == null)
        currentIndicator = handleMissing(currentChunk);
      result[chunkIndex] = currentIndicator;
    }
    return result;
  }
  
  private double[] handleMissing(String currentChunk)
  { 
    BriefLog.warnOnce("WARNING: symbol " + currentChunk + " is unknown (setting to unknown)");
    return getUnknownIndicator();
  }

  public double[] getUnknownIndicator()
  {
    if (_unknownIndicator == null)
    {
      _unknownIndicator = new double[nSymbols()];
      for (int symbolIndex = 0; symbolIndex < nSymbols(); symbolIndex++)
        _unknownIndicator[symbolIndex] = 1.0;
    }
    return _unknownIndicator;
  }
  
  public int nSymbols() 
  {
    return orderedSymbols.size();
  }
  
  public Map<String,double[]> getIndicators()
  {
    if (_indicators == null)
    {
      _indicators = Maps.newHashMap();
      Indexer<String> indexer = getIndexer();
      for (String symbol : orderedSymbols)
      {
        checkChunkLength(symbol.length());
        double [] current = new double[nSymbols()];
        current[indexer.o2i(symbol)] = 1.0;
        BriefCollections.putNoOverwrite(_indicators, symbol.toUpperCase(), current);
      }
      for (String ambiguitySymbol : ambiguousSymbols.keySet())
      {
        checkChunkLength(ambiguitySymbol.length());
        double [] current = new double[nSymbols()];
        for (String possibility : ambiguousSymbols.get(ambiguitySymbol))
          current[indexer.o2i(possibility)] = 1.0;
        BriefCollections.putNoOverwrite(_indicators, ambiguitySymbol.toUpperCase(), current);
      }
    }
    return _indicators;
  }
  
  private void checkChunkLength(int length)
  {
    if (chunkLength == null)
      chunkLength = length;
    if (chunkLength != length)
      throw new RuntimeException();
  }

  public Indexer<String> getIndexer()
  {
    if (_indexer == null) 
      _indexer = new Indexer<String>(orderedSymbols);
    return _indexer;
  }
  
  private PhylogeneticObservationFactory(List<String> orderedSymbols,
      Map<String, Set<String>> ambiguousSymbols, boolean caseSensitive)
  {
    this.orderedSymbols = orderedSymbols;
    this.ambiguousSymbols = ambiguousSymbols;
    this.caseSensitive = caseSensitive;
  }
  
  public static PhylogeneticObservationFactory nucleotidesFactory()
  {
    return fromResource("/conifer/io/dna-iupac-encoding.txt");
  }
  
  public static PhylogeneticObservationFactory fromResource(String resourceURL)
  {
    String jsonString = BriefIO.resourceToString(resourceURL); 
    return new Gson().fromJson(jsonString, PhylogeneticObservationFactory.class);
  }

  public static void main(String [] args)
  {
    System.out.println(EJMLUtils.toString(new SimpleMatrix(nucleotidesFactory().site2CharacterIndicators("asldkfjaslkfasdf"))));
  }
}
