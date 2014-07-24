package conifer.io;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;

import briefj.BriefCollections;
import briefj.BriefIO;
import briefj.BriefLog;
import briefj.Indexer;

import com.google.common.collect.Maps;
import com.google.gson.Gson;

import conifer.ctmc.RateMatrices;
import conifer.ctmc.SimpleRateMatrix;
import conifer.ctmc.expfam.RateMtxNames;


/**
 * A PhylogeneticObservationFactory contains the information required to 
 * parse sequence to turn them into indicator variables sitting at the leaves
 * of a phylogenetic likelihood model.
 * 
 * The main ingredient is an indexer between the integer and strings.
 * e.g. 0 - A, 1 - C, ... Or a more involved example could be codons,
 * e.g. 0 - AAA, 1 - AAC, ...
 * Note that all strings in the indexer should have the same length (the 
 * chunk length).
 * 
 * Note also that some letters read in the strings could correspond to several 
 * symbols in the indexer (these are called ambiguous symbols). For example,
 * for nucleotides, 'S' could be either 'C' or 'G' (representing uncertainty
 * in the read).
 * 
 * Terminology: we call a chunk a string that can be either an actual symbol or 
 * an ambiguous symbol.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class PhylogeneticObservationFactory
{
  /**
   * 
   * @return The PhylogeneticObservationFactory corresponding to the standard iupac encodings.
   */
  public static PhylogeneticObservationFactory nucleotidesFactory()
  {
    if (_nucleotideFactory == null)
      _nucleotideFactory = fromResource("/conifer/io/dna-iupac-encoding.txt");
    return _nucleotideFactory;
  }
  
  public static PhylogeneticObservationFactory proteinFactory()
  {
    if(_proteinFactory == null)
      _proteinFactory = fromResource("/conifer/io/protein-iupac-encoding.txt");
    return _proteinFactory;
   }
  
  public static PhylogeneticObservationFactory proteinPairFactory()
  {
    if(_proteinPairFactory == null)
      _proteinPairFactory = fromResource("/conifer/io/proteinPair-iupac-encoding.txt");
    return _proteinPairFactory;
   }
  
  public static PhylogeneticObservationFactory selectedFactory(final RateMtxNames selectedRateMtx)
  {
    PhylogeneticObservationFactory result = null;
    if (selectedRateMtx == null) {
      throw new IllegalArgumentException("model is null!");
    }
    
    return selectedRateMtx.getFactory();
    
  }
  
  /**
   * Reads the specifications of a PhylogeneticObservationFactory from a JSON file.
   * 
   * @see src/main/resources/conifer/io/dna-iupac-encoding.txt for an example of the format.
   * 
   * @param jsonFile
   * @return
   */
  public static PhylogeneticObservationFactory fromJSONFile(File jsonFile)
  {
    String jsonString = BriefIO.fileToString(jsonFile);
    return fromJSONString(jsonString);
  }
  
  /**
   *
   * @param sequence To be chunked and indexed.
   * @return An array with rows indexing sites and columns indexing actual symbols, filled with one
   *         if the actual symbol is permitted at that site or zero otherwise.
   */
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

  /**
   * 
   * @return The number of symbols permitted in the observations.
   */
  public int nSymbols() 
  {
    return orderedSymbols.size();
  }
  
  private int getChunkLength()
  {
     return orderedSymbols.get(0).length();
  }
  
  
  /**
   * 
   * @return
   */
  
  public Indexer<String> getIndexer()
  {
    if (_indexer == null) 
      _indexer = new Indexer<String>(orderedSymbols);
    return _indexer;
  }
  
  /**
   * 
   * @return A map from chunk to arrays filled with one
   *         if this actual symbol is permitted at that site or zero otherwise.
   */
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
 
  private final List<String> orderedSymbols;
  private final Map<String, Set<String>> ambiguousSymbols;
  private final boolean caseSensitive;
  
  private transient Integer chunkLength = null;
  private transient Indexer<String> _indexer;
  private transient Map<String,double[]> _indicators;
  private transient double[] _unknownIndicator;
  public int nChunks;
  
  private double[] getUnknownIndicator()
  {
    if (_unknownIndicator == null)
    {
      _unknownIndicator = new double[nSymbols()];
      for (int symbolIndex = 0; symbolIndex < nSymbols(); symbolIndex++)
        _unknownIndicator[symbolIndex] = 1.0;
    }
    return _unknownIndicator;
  }
  
  private static PhylogeneticObservationFactory fromResource(String resourceURL)
  {
    String jsonString = BriefIO.resourceToString(resourceURL); 
    return fromJSONString(jsonString);
  }
  
  private static PhylogeneticObservationFactory fromJSONString(String jsonString)
  {
    return new Gson().fromJson(jsonString, PhylogeneticObservationFactory.class);
  }
  
  private double[] handleMissing(String currentChunk)
  { 
    BriefLog.warnOnce("WARNING: symbol " + currentChunk + " is unknown (setting to unknown)");
    return getUnknownIndicator();
  }

  private void checkChunkLength(int length)
  {
    if (chunkLength == null)
      chunkLength = length;
    if (chunkLength != length)
      throw new RuntimeException();
  }
  
  private PhylogeneticObservationFactory(List<String> orderedSymbols,
      Map<String, Set<String>> ambiguousSymbols, boolean caseSensitive)
  {
    this.orderedSymbols = orderedSymbols;
    this.ambiguousSymbols = ambiguousSymbols;
    this.caseSensitive = caseSensitive;
  }
  
  private static PhylogeneticObservationFactory _nucleotideFactory = null;
  private static PhylogeneticObservationFactory _proteinFactory = null;
  private static PhylogeneticObservationFactory _proteinPairFactory=null;
  
  public int nSites()
  {
      return nChunks;
  }

//  public int nSites()
//  {
//    return BriefCollections.pick(getIndicators().values()).length;
//    //return BriefCollections.pick(getIndicators().values()).length()/factory.getChunkLength());
//  }
}
