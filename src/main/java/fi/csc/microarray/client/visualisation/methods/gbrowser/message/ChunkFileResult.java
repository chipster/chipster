package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;

/**
 * Result containing file contents. Used by the data source layer to return content to the processing (tree) layer. 
 *
 */
public class ChunkFileResult {

	public ChunkFileRequest request;
	public ByteRegion exactRegion;
	public FileParser chunkParser;
	public FsfStatus status;
	public Chunk chunk;

	/**
	 * @param fileRequest
	 * @param avg
	 * @param requestQueueSize
	 *            only to update user interface
	 */
	public ChunkFileResult(Chunk chunk, ChunkFileRequest fileRequest, FileParser inputParser, ByteRegion exactRegion, FsfStatus status) {
		this.request = fileRequest;
		this.chunkParser = inputParser;
		this.exactRegion = exactRegion;
		this.status = status;
		this.chunk = chunk;
	}

}