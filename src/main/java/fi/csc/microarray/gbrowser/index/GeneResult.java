package fi.csc.microarray.gbrowser.index;

import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class GeneResult extends AreaResult {
	
	private String searchString;

	public GeneResult(FsfStatus status, List<RegionContent> contents, String searchString) {
		super(status, contents);
		this.searchString = searchString;
	}

	public Object getSearchString() {
		return searchString;
	}

	public Region getGeneLocation() {
		if (super.getContents().size() == 1) {
			return super.getContents().get(0).region;
		}
		return null;
	}
}
