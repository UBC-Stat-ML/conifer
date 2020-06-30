package conifer.io;

import java.io.File;

import java.util.*;

import briefj.BriefCollections;
import briefj.BriefIO;
import briefj.BriefLog;
import briefj.Indexer;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.gson.Gson;
import com.rits.cloning.Immutable;

import blang.inits.DesignatedConstructor;
import blang.inits.Input;
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
@Immutable
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
  
  
  public static Set<String> mapAmbiguousCodonsToAllPossibleCodons(String codons){
	  
	  // check if the input codon consists of three DNA letters, if it is not three digits, throw an exception
	  if(codons.length()!=3)
		  throw new RuntimeException("The length of the codon string is not three");
	  
	  // get the ambiguous symbols from DNA
	  PhylogeneticObservationFactory dnaFactory = nucleotidesFactory();
	  List<String> dnaStates = dnaFactory.orderedSymbols;
	  
	  Map<String, Set<String>> ambiguousDNAMap = dnaFactory.ambiguousSymbols;
	  Set<String> allAmbiguousDNASymbols = ambiguousDNAMap.keySet();
	  // if there is no ambiguous symbols, return itself,
	  // if there exists an ambiguous symbol, return all possible codons, but remove the duplicates
	  List<Set<String>> result = Lists.newArrayList(new HashSet<>());
	  
	  for(int i=0; i < codons.length(); i++){
		  String element = String.valueOf(codons.charAt(i));
		  boolean eleIsAmbiguousSymbol = allAmbiguousDNASymbols.contains(element);
		  if(!dnaStates.contains(element) && !eleIsAmbiguousSymbol)
			  throw new RuntimeException("The character of the codon is not valid");
		  // loop over all the keys of the ambiguousSymbols
		  if(eleIsAmbiguousSymbol){
			  result.add(i, ambiguousDNAMap.get(element));
		  }else{
			  Set<String> storeElement = new HashSet<>();
			  storeElement.add(element);
			  result.add(i, storeElement);
		  }
	   }
	  
	   Set<String> possibleCodons = new HashSet<>();
	   Set<String> possibleCodonsContainer = new HashSet<>();
	   // enumerate all combinations of the possible codons based on the string of the three positions
	   int iter = 0;
	   while(iter < codons.length()){
		   
		   Set<String> allSymbolsCurPosition = result.get(iter);
		   if(possibleCodons.isEmpty()){
			   possibleCodons.addAll(allSymbolsCurPosition);
		   }else{
			   List<String> backupList = Lists.newArrayList();
			   int count = possibleCodons.size();
			   while(backupList.size() < count){
				   for(String element:possibleCodons){
					   for(String elementOfCurPosition: allSymbolsCurPosition){
						   String combinedElement = element.concat(elementOfCurPosition);
						   possibleCodonsContainer.add(combinedElement);					   
					   		}	
					   backupList.add(element);
				   }
			   }
			   possibleCodons.clear();
			   possibleCodons.addAll(possibleCodonsContainer);
			   possibleCodonsContainer.clear();
		   }
		   iter++;
	   }
	   
	   return possibleCodons;
  }
  
  public static PhylogeneticObservationFactory codonFactory(){
	  
	  boolean caseSensitive = false;
	  List<String> compressedCodonCodes = Lists.newArrayList("GCN", "CGN", "MGR", "AAY", "GAY", "TGY", "CAR", "GAR", "GGN", "CAY",
			  "ATH", "YTR", "CTN", "AAR", "TTY", "CCN", "TCN", "AGY", "ACN", "TAY", "GTN");
	  Map<String, Set<String>> ambiguousSymbols = new HashMap<>();
	  for(String element : compressedCodonCodes){
		  ambiguousSymbols.put(element, mapAmbiguousCodonsToAllPossibleCodons(element));
	  }
	  
	  List<String> allDNAStates = nucleotidesFactory().orderedSymbols;
	  List<String> stoppingCodons = Lists.newArrayList("TAA", "TAG", "TGA");
	  
	  List<String> orderedSymbols = Lists.newArrayList();
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
  
  @DesignatedConstructor
  public static PhylogeneticObservationFactory parse(@Input(formatDescription = "DNA, protein, codon or path to JSON spec") String description)
  {
    String cleanedDescr = description.trim().toUpperCase();
    if (cleanedDescr.equals("DNA"))
      return nucleotidesFactory();
    else if (cleanedDescr.equals("PROTEIN"))
      return proteinFactory();
    else if (cleanedDescr.equals("CODON"))
      return codonFactory();
    else
      return fromJSONFile(new File(description));
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
  public int nSites()
  {
    return BriefCollections.pick(getIndicators().values()).length;
    //return BriefCollections.pick(getIndicators().values()).length()/factory.getChunkLength());
  }

  public static void main(String[] args) {

    PhylogeneticObservationFactory factory = PhylogeneticObservationFactory.codonFactory();
    boolean caseSensitive = factory.caseSensitive;
    System.out.println(caseSensitive);
    
    // check if the function is correct when we input a codon that doesn't have ambiguous characters
    System.out.println(mapAmbiguousCodonsToAllPossibleCodons("AGA"));
    

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
