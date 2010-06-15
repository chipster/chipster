package fi.csc.microarray.databeans.features;

import fi.csc.microarray.databeans.Dataset;

public interface FeatureProvider {
	
	public void setName(String name);
	public String getName();
	public Feature createFeature(String namePostfix, Dataset bean);
}
