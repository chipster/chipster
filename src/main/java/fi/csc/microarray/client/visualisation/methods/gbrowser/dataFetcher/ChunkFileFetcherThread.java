package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.IOException;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ByteRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ChunkFileRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ChunkFileResult;

/**
 * <p>Chunk file fetcher is data fetching implementation with a tree type 
 * data structure and it was made primarily for tab separated text files. This class is the 
 * lowest type of file reading. This is done in separate thread to avoid any other things to slow 
 * down file reading.</p> 
 * 
 * <p>Communication between the threads is done with the queues. {@link ChunkFileRequest} are file locations that
 * need reading and content of the file is returned with {@link ChunkFileResult} objects.</p>
 * 
 * @author Petri Klemelä
 */
public class ChunkFileFetcherThread extends Thread {

	private BlockingQueue<ChunkFileRequest> fileRequestQueue;
	private ConcurrentLinkedQueue<ChunkFileResult> fileResultQueue;

	private ChunkTreeHandlerThread areaRequestThread;

	private ChunkDataSource dataSource;

	private FileParser inputParser;
	
	private boolean poison = false;

	public ChunkFileFetcherThread(BlockingQueue<ChunkFileRequest> fileRequestQueue,
	        ConcurrentLinkedQueue<ChunkFileResult> fileResultQueue,
	        ChunkTreeHandlerThread areaRequestThread, FileParser inputParser) {

		this.fileRequestQueue = fileRequestQueue;
		this.fileResultQueue = fileResultQueue;
		this.areaRequestThread = areaRequestThread;
		this.inputParser = inputParser;

		this.setDaemon(true);

		this.dataSource = areaRequestThread.getFile();
	}

	public void run() {

		while (!poison) {
			try {
				processFileRequest(fileRequestQueue.take());
				
			} catch (IOException e) {
				e.printStackTrace(); // FIXME fix exception handling
			} catch (InterruptedException e) {
				e.printStackTrace(); // FIXME fix exception handling
			}
		}
		
		dataSource.close();
	}
	
	/**
	 * Reads the requested parts of the file and returns them with FileResult objects. There are 
	 * two ways of reading, one to use when the location of the line changes isn't known and other
	 * to be used when the same location is read later and the exact location of lines is known
	 * already.
	 * 
	 * This method assumes file content to be separated with new line characters. If the lines
	 * are long (more than 1000 bytes) the current implementation may lose the last line.
	 * 
	 * @See FsfStatus for description of chunk
	 * @param fileRequest
	 * @throws IOException
	 */
	private void processFileRequest(ChunkFileRequest fileRequest) throws IOException {
		
		if (fileRequest.status.poison) {
			poison = true;
			return;
		}
		
		Chunk chunk = new Chunk();
		ByteRegion exactRegion = null;
		
		
		/* If the fileRequest.byteRegion.exact is set, the requested area starts from the beginning 
		 * of the line and ends to the new line character of the same or other line.
		 */
		if (fileRequest.byteRegion.exact) {
			
			// FIXME This is never used

			byte[] byteChunk = new byte[(int)fileRequest.byteRegion.getLength()];
				
			dataSource.read(fileRequest.byteRegion.start, byteChunk);			
			
			chunk.setContent(new String(byteChunk));

		/* fileRequest.byteRegion.exact isn't set and the location of line changes isn't known. 
		 * The returned chunk should contain only full lines starting from the line just after
		 * the first new line character after the request start location. The last returned line
		 * should be the first whose ending new line comes after the end of the request. This way
		 * no lines will be lost between the chunks, even though the exact byte location of the 
		 * new line characters isn't known.   
		 */
		} else {
			
			//FIXME There shouldn't be other limits for String length than Integer.MAX_VALUE and
			//memory heap size, but there seems to be some problems when if the length of chunks is
			//bigger than a couple thousand bytes.
			
			// some extra to get the last line fully
			byte[] byteChunk = new byte[(int)fileRequest.byteRegion.getLength() + 1000];
			
			int length = dataSource.read(fileRequest.byteRegion.start, byteChunk);
			
			String file = new String(byteChunk).substring(0, length);		
			
			exactRegion = new ByteRegion();
			int i = 0;
			
			if(fileRequest.byteRegion.start != 0) {
				for (; ; i++) {
					
					if ( i >= file.length()) {
						//not a single new line found, source file is broken
						return;
					}
					
					if (file.charAt(i) == '\n') {
						i++;
						exactRegion.start = fileRequest.byteRegion.start + i;
						break;
					}
				}
			} else {
				exactRegion.start = 0l;
			}
			
			StringBuffer lines = new StringBuffer();
			
			for (; ; i++) {
				
				
				lines.append(file.charAt(i));		
				
				if (file.charAt(i) == '\n' && i > fileRequest.byteRegion.getLength()) {
					break;
				}
				
				if ( i >= file.length() - 1) {	
					
					// buffer ended before the new line character, discard the last line 
					lines.setLength(lines.lastIndexOf("\n") + 1);
					break;
				}
			}
			
			exactRegion.end = fileRequest.byteRegion.start + i;			
			exactRegion.exact = true;

			chunk.setContent(lines.toString());
			chunk.setByteLocation(exactRegion.start);
		}

		fileRequest.status.maybeClearQueue(fileResultQueue);
		fileRequest.status.fileRequestCount = fileRequestQueue.size();

		ChunkFileResult result = new ChunkFileResult(chunk, fileRequest, inputParser, exactRegion, fileRequest.status);

		fileResultQueue.add(result);
		areaRequestThread.notifyAreaRequestHandler();

	}

	public long getFileLength() {
		if (this.isAlive()) {
			//This requirement isn't really obligatory, but to avoid problems all the communication
			//between threads is done trough the message queues
			throw new IllegalStateException("must be called before the thread is started");
		}

		try {
			return dataSource.length();
		} catch (IOException e) {
			e.printStackTrace(); // FIXME fix exception handling 
		}
		return 0;
	}
	
	public String toString() {
		return this.getClass().getName() + " - " + dataSource;
	}
}
