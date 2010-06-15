package fi.csc.microarray.databeans.features;


import fi.csc.microarray.databeans.Dataset;
import fi.csc.microarray.exception.MicroarrayException;

public abstract class BasicFeature extends FeatureBase implements Feature {
	
	private Dataset bean;
	private FeatureProvider factory;
	
	protected BasicFeature(Dataset bean, FeatureProvider factory) {
		this.bean = bean;
		this.factory = factory;
	}
	
	protected Dataset getDataBean() {
		return bean;
	}
	
	public Iterable<Float> asFloats() throws MicroarrayException {
		return null;
	}
	
	public Table asTable() throws MicroarrayException {
		return null;
	}
	
	public Iterable<String> asStrings() throws MicroarrayException {
		return null;		
	}
	
	public boolean exists() {
		try {
			// try to get it, if fails, it does not exist
			if (asFloats() != null) {
				return true;
			}
			if (asStrings() != null) {
				return true;
			}
			if (asTable() != null) {
				return true;
			}
		} catch (MicroarrayException e) {
			// ignore
		}
		return false;
	}
	
	public String getName() {
		return factory.getName();
	}
}
