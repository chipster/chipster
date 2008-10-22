package fi.csc.microarray.databeans.features;

import fi.csc.microarray.databeans.DataBean;

public interface FeatureProvider {
	
	public void setName(String name);
	public String getName();
	public Feature createFeature(String namePostfix, DataBean bean);
}
