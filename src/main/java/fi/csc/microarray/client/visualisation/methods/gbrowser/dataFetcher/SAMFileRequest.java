package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;

public class SAMFileRequest {

	public AreaRequest areaRequest;
	private FsfStatus status;
	private BpCoord from;
	private BpCoord to;

	public SAMFileRequest(AreaRequest request, BpCoord from, BpCoord to, FsfStatus status) {
		this.areaRequest = request;
		this.status = status;
		this.from = from;
		this.to = to;
	}

	public BpCoord getFrom() {
		return from;
	}

	public BpCoord getTo() {
		return to;
	}

	public FsfStatus getStatus() {
		return status;
	}
}
