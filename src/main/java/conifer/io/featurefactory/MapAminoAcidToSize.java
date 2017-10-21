package conifer.io.featurefactory;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.internal.Lists;

public class MapAminoAcidToSize implements MapSingleStateToFeatures {
	public List<String> stateSpace = null;
	
	public MapAminoAcidToSize(){
		
	}
	
	@Override
	public Map<String, String> getStatesAndFeatures(){
		
		Map<String, List<String>> sizeToAminoAcid = new HashMap<>();
		sizeToAminoAcid.put("Micro", Lists.newArrayList("A", "G", "S", "P", "C", "T", "N", "D"));
		sizeToAminoAcid.put("Big", Lists.newArrayList("R", "E", "Q", "H", "I", "L", "K", "M", "F", "W", "Y", "V"));
		Map<String, String> aminoAcidToSize = new HashMap<>();
		for(String ele:sizeToAminoAcid.keySet()){
			List<String> values = sizeToAminoAcid.get(ele);
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

	
