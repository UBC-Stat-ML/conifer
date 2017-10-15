package conifer.io;

import java.io.File;
import java.util.*;

import briefj.BriefCollections;
import briefj.BriefIO;
import briefj.BriefLog;
import briefj.Indexer;

import com.beust.jcommander.internal.Lists;
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
 * @author Sohrab Salehi (sohrab.salehi@gmail.com)
 */
public class PhylogeneticObservationFactory
{
  /**
   * 
   * @return The PhylogeneticObservationFactory corresponding to the standard iupac encodings.
   */

  public static PhylogeneticObservationFactory proteinFactory()
  {

    boolean caseSensitive = false;
    List<String> orderedSymbols = Lists.newArrayList("A", "R", "N", "D", "C", "Q", "E", "G", "H", "L", "I", "K", "M", "F", "P", "S", "T", "W", "Y", "V");

    Map<String, Set<String>> ambiguousSymbols = new HashMap<>();
    ambiguousSymbols.put("X", new HashSet<String>(orderedSymbols));
    ambiguousSymbols.put("-", new HashSet<String>(orderedSymbols));
    ambiguousSymbols.put("B", new HashSet<String>(Lists.newArrayList("N", "D")));
    ambiguousSymbols.put("Z", new HashSet<>(Lists.newArrayList("Q", "E")));

    return new PhylogeneticObservationFactory(orderedSymbols, ambiguousSymbols, caseSensitive);

  }
  
  
  public static PhylogeneticObservationFactory codonFactory(){
	  
	  boolean caseSensitive = false;
	  List<String> allDNAStates = Lists.newArrayList("A", "C", "G", "T");
	  List<String> stoppingCodons = Lists.newArrayList("TAA", "TAG", "TGA");
	  
	  List<String> orderedSymbols = Lists.newArrayList();
	  Map<String, Set<String>> ambiguousSymbols = new HashMap<>();
	  
	  // write code to generate all codon states
	  for(String ele1:allDNAStates){
		  for(String ele2:allDNAStates){
			  for(String ele3:allDNAStates){
				  String result = ele1.concat(ele2).concat(ele3);
				  orderedSymbols.add(result);
				  }			  
		  }	  
	   }
	  
	  orderedSymbols.removeAll(stoppingCodons);
	  
	  // no inclusion of keys TRA, since it is compressed code for stopping codon
	  List<String> nPrefix = Lists.newArrayList("GC", "CG", "GG", "CT", "CC", "TC", "AC", "GT");
      List<String> yPrefix = Lists.newArrayList("AA", "GA", "TG", "CA", "TT", "AG", "TA");
      List<String> hPrefix = Lists.newArrayList("AT");
      List<String> rPrefix = Lists.newArrayList("CA", "GA", "AA"); // remove TA since TAR is used for stopping codons
      List<String> multiSymbols = Lists.newArrayList("MGR");

      
      Map<String, List<String>> singleCharMaps = new HashMap<>();
      singleCharMaps.put("N", Lists.newArrayList("A", "C", "G", "T"));
      singleCharMaps.put("Y", Lists.newArrayList("C", "T"));
      singleCharMaps.put("R", Lists.newArrayList("A", "G"));
      singleCharMaps.put("H", Lists.newArrayList("A", "C", "T"));
      
      // create all the keys and values for ambiguous symbols
      for(String ele: nPrefix){
    	  String key = ele.concat("N");
    	  // each key corresponds to four values, for example GCN correspond to GCA, GCC, GCG, GCT
    	  Set<String> values = new HashSet<>();
    	  for(String ele2:singleCharMaps.get("N")){
    		  String value = ele.concat(ele2);
    		  values.add(value);
    	  	}
    	  ambiguousSymbols.put(key, values);
    	  }
      
      for(String ele: yPrefix){
    	  String key = ele.concat("Y");
    	  Set<String> values = new HashSet<>();
    	  for(String ele2:singleCharMaps.get("Y")){
    		  String value = ele.concat(ele2);
    		  values.add(value);    		  
    	  }
    	  ambiguousSymbols.put(key, values);	  
      }
      
      for(String ele: rPrefix){
    	  String key = ele.concat("R");
    	  Set<String> values = new HashSet<>();
    	  for(String ele2:singleCharMaps.get("R")){
    		  String value = ele.concat(ele2);
    		  values.add(value);    		  
    	  }
    	  ambiguousSymbols.put(key, values);	  
      }
      
      
      for(String ele: hPrefix){
    	  String key = ele.concat("H");
    	  Set<String> values = new HashSet<>();
    	  for(String ele2:singleCharMaps.get("H")){
    		  String value = ele.concat(ele2);
    		  values.add(value);
    		  }
    	  ambiguousSymbols.put(key, values);
    	  }
      
      ambiguousSymbols.put("MGR", new HashSet<>(Lists.newArrayList("AGA", "AGG", "CGA", "CGG")));
      ambiguousSymbols.put("YTR", new HashSet<>(Lists.newArrayList("CTA", "CTG", "TTA", "TTG")));
      
      return new PhylogeneticObservationFactory(orderedSymbols, ambiguousSymbols, caseSensitive);
  }

  public static PhylogeneticObservationFactory nucleotidesFactory(){

      boolean caseSensitive = false;
      List<String> orderedSymbols = Lists.newArrayList("A", "C", "G", "T");
      Map<String, Set<String>> ambiguousSymbols = new HashMap<>();
      Set<String> wholeValue = new HashSet<String>(Lists.newArrayList("A", "C", "G", "T"));
      Map<String, Set<String>> allMaps = new HashMap<>();

      List<String> eleHasAppeared = Lists.newArrayList();
      // create all possible values for ambiguousSymbols
      for(String element:wholeValue){
        List<String> states = Lists.newArrayList("A", "C", "G", "T");
        states.remove(element);
        List<String> complementSet = states;

        // putting all the three states into allMaps
        Set<String> threeLetterStates = new HashSet<String>(complementSet);
        String threeLetterKeyValues = "";
        for(String secondElement: complementSet){
          threeLetterKeyValues = threeLetterKeyValues.concat(secondElement);
        }

        allMaps.put(threeLetterKeyValues, threeLetterStates);
        
        for(String secondElement : complementSet){
        	Set<String> result = new HashSet<>(Lists.newArrayList(element));
            result.add(secondElement);
            String keyValue = element.concat(secondElement);
            allMaps.put(keyValue, result);
            }  
        	     
      }

    Set<String> tValue = new HashSet<String>(Lists.newArrayList("T"));
    ambiguousSymbols.put("?", wholeValue);
    ambiguousSymbols.put("-", wholeValue);
    ambiguousSymbols.put("N", wholeValue);
    ambiguousSymbols.put("U", tValue);
    ambiguousSymbols.put("W", allMaps.get("AT"));
    ambiguousSymbols.put("S", allMaps.get("CG"));
    ambiguousSymbols.put("M", allMaps.get("AC"));
    ambiguousSymbols.put("K", allMaps.get("GT"));
    ambiguousSymbols.put("R", allMaps.get("AG"));
    ambiguousSymbols.put("Y", allMaps.get("CT"));
    ambiguousSymbols.put("B", allMaps.get("CGT"));
    ambiguousSymbols.put("D", allMaps.get("AGT"));
    ambiguousSymbols.put("H", allMaps.get("ACT"));
    ambiguousSymbols.put("V", allMaps.get("ACG"));

    return new PhylogeneticObservationFactory(orderedSymbols, ambiguousSymbols, caseSensitive);

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
  
  public int getChunkLength()
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
 
  /**
   * 
   * @return Inverse of getIndicators(), with string value of the indicator arrays as keys and the 
   * string chunk as values.
   *  
   */
  public Map<String,String> getIndicator2ChunkMap() 
  {
	  Map<String, String> a2s = Maps.newHashMap();
	  // TODO: this will not include U, as it comes after T, what could be done for not 1-to-1 relations?
	  for (Map.Entry<String, double[]> e : this.getIndicators().entrySet()) {
		  if (e.getKey() != "U")
			  a2s.put(Arrays.toString(e.getValue()), e.getKey());
	  }
	  
	  return a2s;
  }
  
  public final List<String> orderedSymbols;
  public final Map<String, Set<String>> ambiguousSymbols;
  public final boolean caseSensitive;
  
  private transient Integer chunkLength = null;
  private transient Indexer<String> _indexer;
  private transient Map<String,double[]> _indicators;
  private transient double[] _unknownIndicator;
  
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
  
  public PhylogeneticObservationFactory(List<String> orderedSymbols,
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
    return BriefCollections.pick(getIndicators().values()).length;
    //return BriefCollections.pick(getIndicators().values()).length()/factory.getChunkLength());
  }

  public static void main(String[] args) {

    PhylogeneticObservationFactory factory = PhylogeneticObservationFactory.codonFactory();
    boolean caseSensitive = factory.caseSensitive;
    System.out.println(caseSensitive);

    System.out.println(Arrays.deepToString(factory.orderedSymbols.toArray()));
    System.out.println(factory.orderedSymbols.size());

    // print all the keys and values of a HashMap
    for(String name: factory.ambiguousSymbols.keySet()){
      String key = name.toString();
      Set<String> value = factory.ambiguousSymbols.get(key);
      System.out.println(key+" "+ Arrays.deepToString(value.toArray()));
    }
  }

}
