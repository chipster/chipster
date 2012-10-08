package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.IOException;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;

import org.broad.tribble.readers.TabixReader;
import org.broad.tribble.readers.TabixReader.Iterator;

import fi.csc.microarray.client.visualisation.methods.gbrowser.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

/**
 * @author Aleksi Kallio, Petri Klemelä
 *
 */
public abstract class TabixFileFetcherThread extends Thread {

	protected BlockingQueue<BpCoordFileRequest> fileRequestQueue;
	protected ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue;

	protected TabixDataSource dataSource;

	protected boolean poison = false;

	public void run() {

		while (!poison) {

			try {
				
				for (BpCoordFileRequest fileRequest : fileRequestQueue) {
					if (fileRequest.getStatus().poison) {
						poison = true;
					}
				}
				
				processFileRequest(fileRequestQueue.take());

			} catch (IOException e) {
				e.printStackTrace(); // FIXME fix exception handling
			} catch (InterruptedException e) {
				e.printStackTrace(); // FIXME fix exception handling
			}
		}
	}

	abstract protected void processFileRequest(BpCoordFileRequest take) throws IOException;

	protected Iterator getTabixIterator(Region request) {
		String chromosome = request.start.chr.toNormalisedString();

		//limit to integer range
		int start = (int) Math.min(Integer.MAX_VALUE, request.start.bp);
		int end = (int) Math.min(Integer.MAX_VALUE, request.end.bp);
		
		//Extend area to be able to draw introns at screen edge, but don't go over MAX_VALUE, or below 1
		//TODO Be more clever to avoid getting so much useless data
		int EXTRA = 500000; //O,5M should be enought for the longest human introns http://www.bioinfo.de/isb/2004040032/
		
		start = (int) Math.max((long)start - EXTRA, 1);
		end = (int) Math.min((long)end + EXTRA, Integer.MAX_VALUE);

		//Check that region is below max bin size of Tabix
		int MAX_BIN_SIZE = 512*1024*1024 - 2;

		start = (int) Math.min(MAX_BIN_SIZE, start);
		end = (int) Math.min(MAX_BIN_SIZE, end);

		start = (int) Math.max(1, start);
		end = (int) Math.max(1, end);

		String queryRegion = chromosome + ":" + start + "-" + end;

		TabixReader.Iterator iter = null;
		
		try {
			iter = dataSource.getReader().query(queryRegion);
		} catch (ArrayIndexOutOfBoundsException e) {
			//No such chromosome
		}
			
		return iter;
	}

	public String toString() {
		return this.getClass().getName() + " - " + dataSource;
	}
}