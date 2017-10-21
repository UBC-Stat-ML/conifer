package conifer.io.featurefactory;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class MapStateSpaceCharacterToGTR implements MapSingleStateToFeatures {

	public List<String> stateSpace = null;
	public MapStateSpaceCharacterToGTR(){
		
		
	}
	
	@Override
	public Map<String, String> getStatesAndFeatures() {
		// TODO Auto-generated method stub
		Map<String, String>  charToGTR = new LinkedHashMap<>();
		for(String ele: this.stateSpace){
			charToGTR.put(ele, ele);
		}		
		return charToGTR;
	}

	@Override
	public void setStateSpace(List<String> stateSpace) {
		// TODO Auto-generated method stub
		this.stateSpace = stateSpace;
	}

}
