package conifer.io;

import static blang.variables.RealVariable.real;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;

import blang.annotations.FactorArgument;
import blang.annotations.Processors;
import blang.annotations.Samplers;
import blang.variables.RealVariable;
import briefj.Indexer;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import conifer.Parsimony;
import conifer.ParsimonyVector;
import conifer.TreeNode;
import conifer.models.CNPair;
import conifer.models.CNSpecies;
import conifer.moves.CopyNumberTreeSampler;
import conifer.processors.CopyNumberTreeProcessor;

@Samplers({ CopyNumberTreeSampler.class })
@Processors({ CopyNumberTreeProcessor.class })
/**
 * Keeps the raw count data as CNSpacies,
 * 
 * @author Sohrab Salehi (sohrabsalehi@gmail.com)
 * @author Sean Jewell (jewellsean@gmail.com)
 */
public class CopyNumberTreeObservation implements TreeObservations {
	
    @FactorArgument
    public final RealVariable betaBinomialprecision = real(1);
    
    public final Parsimony parsimony;
    
    private Map<String, Integer> leafOrder = null;
    private Map<Integer, String> leafString = null;
    private Map<String, Set<CNSpecies>> cnSpeciesClusters = new LinkedHashMap<String, Set<CNSpecies>>();
	
    // raw count data E(v)
    private final Set<CNSpecies> cnSpecies = new LinkedHashSet<CNSpecies>();

    // Y(v): v \in L
    private LinkedHashMap<TreeNode, double[][]> currentCTMCState = Maps.newLinkedHashMap();

    private final int nSites;
    
    private final Set<TreeNode> leaves;
    
    public Map<String, Integer> getLeafOrder()
	{
	    if (leafOrder != null)
	    {
	        return leafOrder;
	    }
	    
	    Map<String, Set<CNPair>> emissions = getEmissionAtSite(0);
	    leafOrder = new HashMap<String, Integer>();
	    Integer i = 0; 
	    for (String s : emissions.keySet())
	        leafOrder.put(s, i++);
	    return leafOrder; 
	}
	
	private Map<Integer, String> getLeafName()
	{
	    if (leafString != null)
	        return leafString;
	    
	    Map<Integer, String> leafString = new HashMap<Integer, String>();
	    for(String s : leafOrder.keySet())
	        leafString.put(leafOrder.get(s), s);

	    return leafString;
	}
	
	public String getLeafString(Integer leaf)
	{
	    return getLeafName().get(leaf);
	}
	
	public int getLeafOrder(String s)
	{
	    return getLeafOrder().get(s).intValue();
	}
	

	public CopyNumberTreeObservation(Set<CNSpecies> cnSpecies) {
		nSites = cnSpecies.iterator().next().getCnPairs().size();
		setCNSpecies(cnSpecies);
		LinkedHashMap<TreeNode, List<CNPair>> leaves = CNParser.getNodeMap(cnSpecies);
		leaves.put(TreeNode.withLabel("root"), null);
		this.leaves = leaves.keySet();
		populateClusterEmissions();
		this.parsimony = initalizeParsimony();
		initalizeLeafStates();    
	}

	public CopyNumberTreeObservation(int nSites, Set<TreeNode> leaves)
    {
	    this.nSites = nSites; 
	    this.leaves = leaves; 
	    this.parsimony = initalizeParsimony();
	    initalizeLeafStates();
    }

    public Set<TreeNode> getLeaves()
	{
	    return leaves;
	}
	
	private Parsimony initalizeParsimony()
	{
	    this.leafOrder = getLeafOrder();
	    this.leafString = getLeafName();
	    return new Parsimony(ParsimonyVector.oneInit(nSites, leafString));
	}
	
	private void initalizeLeafStates()
	{
	    double uniform[][] = new double[nSites][Indexers.copyNumberCTMCIndexer().size()];
	    for (int i = 0; i < nSites; i++)
	        bayonet.opt.DoubleArrays.initialize(uniform[i], 1);
	    for (TreeNode t : getLeaves())
	    {
	        if (t.toString() == "root")
	            currentCTMCState.put(t, rootLeaf());    
	        currentCTMCState.put(t, uniform);
	    }
	}

	private double[][] rootLeaf()
	{
	    double[] rootLeaf = new double[Indexers.copyNumberCTMCIndexer().size()];
	    double[][] root = new double[nSites][Indexers.copyNumberCTMCIndexer().size()];
	    rootLeaf[2] = 1;
	    for (int i = 0; i < nSites; i++)
	        root[i] = rootLeaf;
	    return root;
	}

	@Override
	public List<TreeNode> getObservedTreeNodes() {
		List<TreeNode> result = Lists.newArrayList();
		for (TreeNode key : currentCTMCState.keySet())
			result.add(key);
		return result;
	}

	@Override
	public Object get(TreeNode leaf) {
		return currentCTMCState.get(leaf);
	}

	public double[] getSite(TreeNode leaf, int site) {
		double[][] ctmcStateSpace = currentCTMCState.get(leaf);
		return ctmcStateSpace[site];
	}

	public void setSite(TreeNode leaf, int site, double[] ctmcUpdate) {
		double[][] ctmc = currentCTMCState.get(leaf);
		ctmc[site] = ctmcUpdate;
		currentCTMCState.put(leaf, ctmc);
	}

	@Override
	public void set(TreeNode leaf, Object data) {
		double[][] cast = (double[][]) data;
		if (cast.length != nSites)
			throw new RuntimeException();
		currentCTMCState.put(leaf, cast);
	}

	@Override
	public void clear() {
		getCnSpecies().clear();
		currentCTMCState.clear();
	}

	public void clearCurrentState() {
		currentCTMCState.clear();
	}

	@Override
	public int nSites() {
		return nSites;
	}

	public Map<TreeNode, List<CNPair>> getTreeNodeRepresentation() {
		return CNParser.getNodeMap(getCnSpecies());
	}

	/**
	 * @return the cnSpecies
	 */
	public Set<CNSpecies> getCnSpecies() {
		return cnSpecies;
	}

	public void setCNSpecies(Set<CNSpecies> cnSpecies) {
		for (CNSpecies s : cnSpecies) {
			this.cnSpecies.add(s);
		}
	}

	public Map<String, Set<CNPair>> getEmissionAtSite(int i)
	{
	    Map<String, Set<CNPair>> emissions = new HashMap<String, Set<CNPair>>();

	    for (String cluster : cnSpeciesClusters.keySet())
	    {
	        Set<CNSpecies> clusterSpecies = cnSpeciesClusters.get(cluster);
	        Set<CNPair> cnClusterPairs = new HashSet<CNPair>();
	        Iterator<CNSpecies> itPairs = clusterSpecies.iterator();

	        while(itPairs.hasNext())
	        {
	            CNSpecies pairs = itPairs.next();
	            cnClusterPairs.add(pairs.getCnPairs().get(i));
	        }
	        emissions.put(cluster, cnClusterPairs);
	    }

	    return emissions; 
	}
	
	public List<String> getClusters(Set<CNSpecies> data)
	{
	    List<String> clusters = new ArrayList<String>();
	    for (CNSpecies cn : data)
	        {
	            if (!clusters.contains(cn.getClusterID()))
	                clusters.add(cn.getClusterID());
	        }
	    return clusters; 
	}
	
	private void populateClusterEmissions()
	{
	    for (CNSpecies s : cnSpecies)
	    {
	        String cluster = s.getClusterID();
	        Set<CNSpecies> clusterElements;
	        if (cnSpeciesClusters.containsKey(cluster))
	            clusterElements = cnSpeciesClusters.get(cluster);
	        else
	            clusterElements = new HashSet<CNSpecies>();
	        clusterElements.add(s);
	        cnSpeciesClusters.put(cluster, clusterElements);
	    }
	    
	}
	

	// Return the corresponding state
	public String CNPairFromIndicatorArray(double[] indicatorArray) {
		Indexer<CNPair> indexer = Indexers.CNPairIndexer();
		return indexer.i2o(ArrayUtils.indexOf(indicatorArray, 1)) + "";
	}

	public String CNListFromCTMCStateSpace(double[][] ctmcStateSpace) {
		StringBuilder result = new StringBuilder();
		for (int i = 0; i < nSites; i++) {
			result.append(CNPairFromIndicatorArray(ctmcStateSpace[i]) + "  ");
		}
		return result.toString();
	}

	@Override
	public String toString() {
		StringBuilder result = new StringBuilder();

		LinkedHashMap<TreeNode, List<CNPair>> emissions = CNParser.getNodeMap(cnSpecies);

		for (TreeNode node : leaves) {
			// CTMC state (copy numbers Y(v): v \in L)
			// result.append(node.toString() + "_[CTMC_STATE]: " +
			// Arrays.deepToString(currentCTMCState.get(node)) + " \n");
			result.append(node.toString() + "_[CTMC_STATE]: " + CNListFromCTMCStateSpace(currentCTMCState.get(node))
					+ " \n");

			// Emission State (count data E(v))
			if (emissions.size() > 0)
				result.append(node.toString() + "_[EMISSION]: " + emissions.get(node).toString() + " ");
		}

		return result.toString();
	}

	// return a clone
	public LinkedHashMap<TreeNode, double[][]> getCurrentCTMCState() {
		return new LinkedHashMap<TreeNode, double[][]>(currentCTMCState);
	}

	public List<CNPair> getCNListFromIndicatorMatrix(double[][] indicatorMatrix) {
	  // (0, 2)
	  List<CNPair> result = Lists.newArrayList();
	  
	  for (int i = 0; i < nSites; i++) {
		  String tempPair = CNPairFromIndicatorArray(indicatorMatrix[i]);
		  tempPair = tempPair.replaceAll("\\s+","");
		  tempPair = tempPair.replace("(", "");
		  tempPair = tempPair.replace(")", "");
		  String[] parts = tempPair.split(",");
		  result.add(new CNPair(Integer.parseInt(parts[0]), 
				  Integer.parseInt(parts[1])));
	  }
	  
	  return result;
  }

	// generate a print-friendly object from currentCTMCState
	public LinkedHashMap<TreeNode, List<CNPair>> getPrintFriendlyCTMCState() {
		LinkedHashMap<TreeNode, List<CNPair>> result = Maps.newLinkedHashMap();

		for (TreeNode node : leaves) {
			result.put(node, getCNListFromIndicatorMatrix(currentCTMCState.get(node)));
		}

		return result;
	}
}
