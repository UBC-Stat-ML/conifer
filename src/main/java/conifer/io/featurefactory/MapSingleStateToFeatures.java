package conifer.io.featurefactory;

import java.util.List;
import java.util.Map;

public interface MapSingleStateToFeatures {
	
	public List<String> stateSpace = null;
	Map<String, String> getStatesAndFeatures();
	void setStateSpace(List<String> stateSpace);
}

