package fi.csc.microarray.databeans.features;

import fi.csc.microarray.exception.MicroarrayException;

public class QueryResult {
	
	private Feature feature;
	
	public QueryResult(Feature feature) {
		this.feature = feature;
	}
	
	public String getName() {
		return feature.getName(); // TODO useless?
	}
	
	public boolean exists() {
		return feature.exists(); // TODO convert to asBoolean 
	}
	
	public Iterable<Float> asFloats() throws MicroarrayException {
		return feature.asFloats();
	}
	
	public Table asTable() throws MicroarrayException {
		return feature.asTable();
	}
	
	public Iterable<String> asStrings() throws MicroarrayException {
		return feature.asStrings();
	}
	
	public String asString() throws MicroarrayException {
		return feature.asString();
	}

	public Float asFloat() throws MicroarrayException {
		return feature.asFloat();
	}

	public Feature asFeature() {
		return feature; // TODO remove, only for code migration
	}
}
