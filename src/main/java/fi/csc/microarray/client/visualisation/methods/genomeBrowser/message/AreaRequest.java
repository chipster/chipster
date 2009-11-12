package fi.csc.microarray.client.visualisation.methods.genomeBrowser.message;

import java.util.Collection;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.Content;


public class AreaRequest extends Region{

	public FsfStatus status;
	public Collection<Content> requestedContents;
	
	public AreaRequest(Region region, Collection<Content> requestedContents, FsfStatus status) {
		super(region.start, region.end);
		this.requestedContents = requestedContents;
		this.status = status;
	}
}
