package conifer.models;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Contains CNPairs and Cluster_ID
 * @author Sohrab Salehi (sohrab.salehi@gmail.com)
 *
 */
// TODO: strongly enforce ordering 
public class CNSpecies {
	
	private String speciesName;
	private String clusterID;
	private List<CNPair> cnPairs;
	
	public CNSpecies(List<CNPair> cnPairs, String clusterID, String speciesName) {
		this.cnPairs = cnPairs;
		this.clusterID = clusterID;
		this.speciesName = speciesName;
	}

	
	/**
	 * @return the clusterID
	 */
	public String getClusterID() {
		return clusterID;
	}

	/**
	 * @return the cnPairs
	 */
	public List<CNPair> getCnPairs() {
		return cnPairs;
	}

	/**
	 * @return the speciesName
	 */
	public String getSpeciesName() {
		return speciesName;
	}
	
	// TODO: maybe add site_id?
	@Override
	public String toString() {
		return "SpeciesName: " + this.speciesName + "\n" + 
	"ClusterID: " + this.clusterID + "\n" + 
				Arrays.deepToString(this.getCnPairs().toArray()) + "\n";
	}

}
