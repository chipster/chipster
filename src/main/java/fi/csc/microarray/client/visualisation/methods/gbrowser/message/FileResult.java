package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;

/**
 * Result containing file contents. Used by the file layer to return content to the tree layer. 
 *
 */
public class FileResult {

	public FileRequest request;
	public ByteRegion exactRegion;
	public FileParser chunkParser;
	public FsfStatus status;
	
	/**
	 * DOCME
	 * 
	 * Chunk is a part of file containing only full lines. The default length of the chunk can be 
	 * defined in the parser classes, but the exact length will be only in the some order of 
	 * magnitude and the exact length depends on how the lines happen to become divided between 
	 * these chunks. 
	 */
	public String chunk;

	/**
	 * @param fileRequest
	 * @param avg
	 * @param requestQueueSize
	 *            only to update user interface
	 */
	public FileResult(String chunk, FileRequest fileRequest, FileParser inputParser, ByteRegion exactRegion, FsfStatus status) {
		this.request = fileRequest;
		this.chunkParser = inputParser;
		this.exactRegion = exactRegion;
		this.status = status;
		this.chunk = chunk;
	}

}