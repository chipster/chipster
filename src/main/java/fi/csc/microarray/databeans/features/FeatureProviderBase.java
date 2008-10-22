package fi.csc.microarray.databeans.features;

public abstract class FeatureProviderBase implements FeatureProvider {

	private String name = null;
	
	public void setName(String name) {
		this.name = name;
	}
	
	public String getName() {
		return name;
	}
}
