package fi.csc.microarray.client.visualisation.methods.genomeBrowser.message;

import java.util.Collection;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ColumnType;


public class AreaRequest extends BpCoordRegion {

	public FsfStatus status;
	public Collection<ColumnType> requestedContents;
	
	public AreaRequest(BpCoordRegion region, Collection<ColumnType> requestedContents, FsfStatus status) {
		super(region.start, region.end);
		this.requestedContents = requestedContents;
		this.status = status;
	}
}
