package conifer.io;

import static blang.variables.RealVariable.real;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import blang.annotations.FactorArgument;
import blang.annotations.Processors;
import blang.annotations.Samplers;
import blang.variables.RealVariable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import conifer.Parsimony;
import conifer.ParsimonyVector;
import conifer.TreeNode;
import conifer.models.CNPair;
import conifer.models.CNSpecies;
import conifer.moves.CopyNumberTreeSampler;
import conifer.processors.CopyNumberTreeProcessor;

@Samplers({CopyNumberTreeSampler.class})
@Processors({CopyNumberTreeProcessor.class})
/**
 * Keeps the raw count data as CNSpacies,
 * 
 * @author Sohrab Salehi (sohrab.salehi@gmail.com)
 * @author Sean Jewell (jewellsean@gmail.com)
 */
public class CopyNumberTreeObservation implements TreeObservations {
	
    @FactorArgument
    public final RealVariable betaBinomialprecision = real(1);
    
    public final Parsimony parsimony;
    
    private Map<String, Integer> leafOrder = null;
	
	public Map<String, Integer> getLeafOrder()
	{
	    if (leafOrder != null)
	    {
	        return leafOrder;
	    }
	    
	    Map<String, CNPair> emissions = getEmissionAtSite(0);
	    leafOrder = new HashMap<String, Integer>();
	    Integer i = 0; 
	    for (String s : emissions.keySet())
	        leafOrder.put(s, i++);
	    return leafOrder; 
	}
	
	public int getLeafOrder(String s)
	{
	    return getLeafOrder().get(s).intValue();
	}
	
	// raw count data E(v)
	private final Set<CNSpecies> cnSpecies = new LinkedHashSet<CNSpecies>();

	// Y(v): v \in L
	private LinkedHashMap<TreeNode, double[][]> currentCTMCState = Maps.newLinkedHashMap();

	private final int nSites;
	
	private final Set<TreeNode> leaves;
	
	public CopyNumberTreeObservation(Set<CNSpecies> cnSpecies) 
	{
		nSites = cnSpecies.iterator().next().getCnPairs().size();
		setCNSpecies(cnSpecies);
		this.parsimony = new Parsimony(ParsimonyVector.oneInit(nSites));
		this.leaves = CNParser.getNodeMap(cnSpecies).keySet();
		initalizeLeafStates();
	}

	public Set<TreeNode> getLeaves()
	{
	    return leaves;
	}
	
	private void initalizeLeafStates()
    {
	    double uniform[][] = new double[nSites][Indexers.copyNumberCTMCIndexer().size()];
	    for (int i = 0; i < nSites; i++)
	        bayonet.opt.DoubleArrays.initialize(uniform[i], 1);
	    for (TreeNode t : getLeaves())
	        currentCTMCState.put(t, uniform);
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
	
	public void setSite(TreeNode leaf, int site, double[] ctmcUpdate)
	{
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

	public Map<String, CNPair> getEmissionAtSite(int i)
    {
        Map<String, CNPair> emissions = new HashMap<String, CNPair>();
        Iterator<CNSpecies> itPairs = cnSpecies.iterator();
        
        while(itPairs.hasNext())
        {
            CNSpecies pairs = itPairs.next();
            emissions.put(pairs.getSpeciesName(), pairs.getCnPairs().get(i));
        }
        
        return emissions; 
    }
}
