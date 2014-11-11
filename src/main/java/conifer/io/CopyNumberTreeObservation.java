package conifer.io;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import conifer.TreeNode;
import conifer.models.CNPair;
import conifer.models.CNSpecies;

/**
 * Keeps the raw count data as CNSpacies,
 * 
 * @author Sohrab Salehi (sohrab.salehi@gmail.com)
 * @author Sean Jewell (jewellsean@gmail.com)
 */
public class CopyNumberTreeObservation implements TreeObservations {
	
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

	public CopyNumberTreeObservation(int nSites) {
		this.nSites = nSites;
	}

	public CopyNumberTreeObservation(Set<CNSpecies> cnSpecies) {
		nSites = cnSpecies.iterator().next().getCnPairs().size();
		setCNSpecies(cnSpecies);
//		initialize();
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
		// TODO: what is the expected behavior, clear all data or just clear the
		// current state
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
	
	
	
	// TODO: complete implementation
//	public void initialize() {
//
//		int nSpecies = cnSpecies.size();
//
//		double[][] data = new double[nSpecies][nSites];
//
//		// initialize the probabilityVector
//		double[] alphas = new double[nCTMCStates];
//		for (int i = 0; i < alphas.length; i++) {
//			alphas[i] = 1;
//		}
//		
//		// initialize currentCTMCState
//		for (CNSpecies s : cnSpecies) {
//			for (int j = 0; j < nSites; j++) {
//				data[j] = Multinomial.generate(new Random(), 1, probabilityVector);
//			}
//			currentCTMCState.put(TreeNode.withLabel(s.getSpeciesName()), data);
//		}

		// TODO: do we have to manually make sure that in the initial state
		// there's at least one state with b=1?...yes, but through aux variable
//	}


	
	// TODO: finish implementation
//	@Override
//	public String toString() {
//
//		StringBuilder result = new StringBuilder();
//		
//		PhylogeneticObservationFactory factory = PhylogeneticObservationFactory.copyNumberCTMCFactory();
//		//Map<String,String> map = factory.getIndicator2ChunkMap();
//		TreeNode t = currentCTMCState.keySet().iterator().next();
//		double[][] d = currentCTMCState.get(t);
////		System.out.println(map.get(Arrays.toString(d[1])));
//		System.out.println(Arrays.toString(d[1]));
//		
//	    //for (TreeNode node : currentCTMCState.keySet())
//	     // result.append(node.toString() + " : " + Arrays.deepToString(currentCTMCState.get(node)) + "\n");
//		return Arrays.toString(d[1]);
////		return "Fixed Data: " + Arrays.deepToString(this.cnSpecies.toArray()) + "\n" + 
////				"Current CTMCStates at leaves: " + result.toString();
//	}
}
