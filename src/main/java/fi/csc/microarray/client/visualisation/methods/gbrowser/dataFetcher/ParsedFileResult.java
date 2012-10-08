package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Result generated from {@link SAMHandlerThread} by {@link SAMFileFetcherThread}.
 * 
 * @author Aleksi Kallio
 *
 */
public class ParsedFileResult {

	private List<RegionContent> content;
	private BpCoordFileRequest fileRequest;
	public AreaRequest areaRequest;
	private FsfStatus status;

	public ParsedFileResult(List<RegionContent> content, BpCoordFileRequest fileRequest, AreaRequest areaRequest, FsfStatus status) {
		this.content = content;
		this.fileRequest = fileRequest;
		this.areaRequest = areaRequest;
		this.status = status;
	}

	public BpCoordFileRequest getFileRequest() {
		return fileRequest;
	}

	public AreaRequest getAreaRequest() {
		return areaRequest;
	}

	public List<RegionContent> getContents() {
		return content;
	}

	public FsfStatus getStatus() {
		return status;
	}

}
