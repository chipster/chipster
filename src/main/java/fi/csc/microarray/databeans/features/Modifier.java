package fi.csc.microarray.databeans.features;

import java.util.List;

public interface Modifier {

	public void setInputs(List<Feature> inputs);
	
	public Feature getOutput();
		
}
