package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.List;


public class GeneResult extends DataResult {
	
	private String searchString;

	public GeneResult(DataStatus status, List<RegionContent> contents, String searchString) {
		super(status, contents);
		this.searchString = searchString;
	}

	public String getSearchString() {
		return searchString;
	}

	public Region getGeneLocation() {
		if (super.getContents().size() == 1) {
			return super.getContents().get(0).region;
		}
		return null;
	}
}
