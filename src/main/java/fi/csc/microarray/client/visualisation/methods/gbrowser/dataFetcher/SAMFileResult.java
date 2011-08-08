package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class SAMFileResult {

	private List<RegionContent> content;
	private SAMFileRequest fileRequest;
	public AreaRequest areaRequest;
	private FsfStatus status;

	public SAMFileResult(List<RegionContent> content, SAMFileRequest fileRequest, AreaRequest areaRequest, FsfStatus status) {
		this.content = content;
		this.fileRequest = fileRequest;
		this.areaRequest = areaRequest;
		this.status = status;
	}

	public SAMFileRequest getFileRequest() {
		return fileRequest;
	}

	public AreaRequest getAreaRequest() {
		return areaRequest;
	}

	public List<RegionContent> getContent() {
		return content;
	}

	public FsfStatus getStatus() {
		return status;
	}

}
