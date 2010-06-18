package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.Collection;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;

public class AreaRequest extends BpCoordRegion {

	private static int MAX_RECURSION_DEPTH = 8;
	
	public FsfStatus status;
	public Collection<ColumnType> requestedContents;
	public int depthToGo;

	public AreaRequest(BpCoordRegion region, Collection<ColumnType> requestedContents, FsfStatus status) {
		this(region, requestedContents, status, MAX_RECURSION_DEPTH);
	}
	public AreaRequest(BpCoordRegion region, Collection<ColumnType> requestedContents, FsfStatus status, int depthToGo) {
		super(region.start, region.end);
		this.requestedContents = requestedContents;
		this.status = status;
		this.depthToGo = depthToGo;
	}
	
	@Override
	public AreaRequest clone() throws CloneNotSupportedException {
		return new AreaRequest(this, this.requestedContents, this.status.clone(), this.depthToGo);
	}
}
