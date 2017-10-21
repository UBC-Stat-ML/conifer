package conifer.io.featurefactory;

import java.util.Map;

import com.google.common.collect.Lists;

import java.util.HashMap;
import java.util.List;

public class MapAminoAcidToPolarity implements MapSingleStateToFeatures {
	
	public List<String> stateSpace = null;
	
	public MapAminoAcidToPolarity(){
		
	}
	
	@Override
	public Map<String, String> getStatesAndFeatures(){
		// Acidic stands for Acidic Polar, Basic stands for Basic Polar, Non stands for NonPolar, Polar stands for Polar
		Map<String, List<String>> polarityToAminoAcid = new HashMap<>();
		polarityToAminoAcid.put("Acidic", Lists.newArrayList("D", "E"));
		polarityToAminoAcid.put("Basic", Lists.newArrayList("R", "H", "K"));
		polarityToAminoAcid.put("Non", Lists.newArrayList("A", "C", "G", "I", "L", "M", "F", "P", "W", "V"));
		polarityToAminoAcid.put("Polar", Lists.newArrayList("N", "Q", "S", "T", "Y"));
		Map<String, String> aminoAcidToPolarity = new HashMap<>();
		for(String ele: polarityToAminoAcid.keySet()){
			List<String> values = polarityToAminoAcid.get(ele);
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
	