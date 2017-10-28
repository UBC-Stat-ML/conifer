package conifer.io.featurefactory;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import com.google.gson.Gson;

import briefj.BriefIO;

public class MapAminoAcidToSize implements MapSingleStateToFeatures {
	public List<String> stateSpace = null;
	
	public MapAminoAcidToSize(){
		
	}
	
	public static Map<String, Set<String>> sizeToAminoAcid = new LinkedHashMap();
	 
	private static Map<String, Set<String>> _sizeToAminoAcid = null;
	
	private static Map<String, Set<String>> fromResource(String resourceURL)
	{
	    String jsonString = BriefIO.resourceToString(resourceURL); 
	    return fromJSONString(jsonString);
	    }
	  
	private static Map<String, Set<String>> fromJSONString(String jsonString)
	{
	    return new Gson().fromJson(jsonString, sizeToAminoAcid.getClass());
	    }
	
	public static Map<String, Set<String>> mapSizeToAminoAcid(){
		
		if (_sizeToAminoAcid == null)
			_sizeToAminoAcid = fromResource("/conifer/io/sizeToAminoAcid.txt");
		    return _sizeToAminoAcid;
		
	}
	
	@Override
	public Map<String, String> getStatesAndFeatures(){
		
		sizeToAminoAcid = mapSizeToAminoAcid();
		Map<String, String> aminoAcidToSize = new HashMap<>();
		for(String ele:sizeToAminoAcid.keySet()){
			Set<String> values = sizeToAminoAcid.get(ele);
			for(String value:values){
				aminoAcidToSize.put(value, ele);
			}
		}
	
		return aminoAcidToSize;
	}
	
	@Override
	public void setStateSpace(List<String> stateSpace) {
		// TODO Auto-generated method stub
		this.stateSpace = stateSpace;
	}
	

}

	
