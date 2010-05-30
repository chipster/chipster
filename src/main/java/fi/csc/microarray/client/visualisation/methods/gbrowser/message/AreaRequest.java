package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.Collection;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;

public class AreaRequest extends BpCoordRegion {

	private static int MAX_RECURSION_DEPTH = 2;
	
	public FsfStatus status;
	public Collection<ColumnType> requestedContents;
	public int depthToGo;

	public AreaRequest(BpCoordRegion region, Collection<ColumnType> requestedContents, FsfStatus status) {
		super(region.start, region.end);
		this.requestedContents = requestedContents;
		this.status = status;
		this.depthToGo = MAX_RECURSION_DEPTH;
	}
	
	@Override
	public AreaRequest clone() throws CloneNotSupportedException {
		return new AreaRequest(this, this.requestedContents, this.status.clone());
	}
}
