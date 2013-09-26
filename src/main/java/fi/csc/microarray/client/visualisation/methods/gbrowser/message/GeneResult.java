package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.List;


public class GeneResult extends DataResult {
	
	private String searchString;

	public GeneResult(DataStatus status, List<Feature> contents, String searchString) {
		super(status, contents);
		this.searchString = searchString;
	}

	public String getSearchString() {
		return searchString;
	}

	public Region getGeneLocation() {
		if (super.getFeatures().size() == 1) {
			return super.getFeatures().get(0).region;
		}
		return null;
	}
}
