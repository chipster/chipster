package fi.csc.microarray.client.visualisation.methods.genomeBrowser.message;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher.ByteChunk;


	
	public class FileResult{
		/**
		 * @param fileRequest
		 * @param avg
		 * @param requestQueueSize only to update user interface
		 */
		public FileResult(FileRequest fileRequest, ByteChunk chunk, FsfStatus status) {
			this.request = fileRequest;
			this.chunk = chunk;
			this.status = status;
		}
		public FileRequest request;
		public ByteChunk chunk;
		public FsfStatus status;
	}