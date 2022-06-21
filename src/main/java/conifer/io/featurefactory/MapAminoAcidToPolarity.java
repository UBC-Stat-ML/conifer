package conifer.io.featurefactory;

import java.util.Map;
import java.util.Set;

import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;

import briefj.BriefIO;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

public class MapAminoAcidToPolarity implements MapSingleStateToFeatures {
	
	public List<String> stateSpace = null;
	
	public MapAminoAcidToPolarity(){
		
	}
	
	public static Map<String, Set<String>> polarityToAminoAcid = new LinkedHashMap();
	 
	private static Map<String, Set<String>> _polarityToAminoAcid = null;
	
	private static Map<String, Set<String>> fromResource(String resourceURL)
	{
	    String jsonString = BriefIO.resourceToString(resourceURL, MapAminoAcidToPolarity.class); 
	    return fromJSONString(jsonString);
	    }
	  
	private static Map<String, Set<String>> fromJSONString(String jsonString)
	{
	    return new Gson().fromJson(jsonString, new TypeToken<Map<String, Set<String>>>(){}.getType());
	    }
	
	public static Map<String, Set<String>> mapPolarityToAminoAcid(){
		
		if (_polarityToAminoAcid == null)
			_polarityToAminoAcid = fromResource("/conifer/io/polarityToAminoAcid.txt");
		    return _polarityToAminoAcid;
		
	}
	
	
	@Override
	public Map<String, String> getStatesAndFeatures(){
		// Acidic stands for Acidic Polar, Basic stands for Basic Polar, Non stands for NonPolar, Polar stands for Polar

		polarityToAminoAcid = mapPolarityToAminoAcid();
		Map<String, String> aminoAcidToPolarity = new HashMap<>();
		for(String ele: polarityToAminoAcid.keySet()){
			Set<String> values = polarityToAminoAcid.get(ele);
			for(String value:values){
				aminoAcidToPolarity.put(value, ele);
			}
					
		}		
		return aminoAcidToPolarity;
	}

	
	@Override
	public void setStateSpace(List<String> stateSpace) {
		// TODO Auto-generated method stub
		this.stateSpace = stateSpace;
	}
	
}	
	