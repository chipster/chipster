package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;

public class SAMFileRequest {

	public AreaRequest areaRequest;
	private FsfStatus status;

	public SAMFileRequest(AreaRequest areaRequest, FsfStatus status) {
		this.areaRequest = areaRequest;
		this.status = status;
	}

	public FsfStatus getStatus() {
		return status;
	}
}
