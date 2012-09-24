package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;

import org.broad.tribble.readers.TabixReader;

import fi.csc.microarray.client.visualisation.methods.gbrowser.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * @author Petri Klemelä
 *
 */
public class BedTabixFileFetcherThread extends TabixFileFetcherThread {

	private BedTabixHandlerThread areaRequestThread;

	public BedTabixFileFetcherThread(
			BlockingQueue<BpCoordFileRequest> fileRequestQueue, 
			ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue, 
			BedTabixHandlerThread areaRequestThread, TabixDataSource dataSource) {

		this.fileRequestQueue = fileRequestQueue;
		this.fileResultQueue = fileResultQueue;
		this.areaRequestThread = areaRequestThread;
		this.dataSource = dataSource;
		this.setDaemon(true);
	}

	@Override
	protected void processFileRequest(BpCoordFileRequest fileRequest) throws IOException {

		if (fileRequest.getStatus().poison) {
			poison = true;
			return;
		}

		List<RegionContent> resultList = new LinkedList<RegionContent>();

		Region region = new Region(fileRequest.getFrom(), fileRequest.getTo());
	
		// Read the given region
		TabixReader.Iterator iter = getTabixIterator(region);

		String line;

		if (iter != null) { //null if there isn't such chromosome in annotations

			while ((line = iter.next()) != null) {
				resultList.add(new RegionContent(parseBedLine(line)));
			}
		}

		ParsedFileResult result = new ParsedFileResult(resultList, fileRequest, fileRequest.areaRequest, fileRequest.getStatus());

		fileResultQueue.add(result);
		areaRequestThread.notifyAreaRequestHandler();	
	}


	private Region parseBedLine(String line) {

		String[] cols;

		cols = line.split("\t");

		String chr = cols[0];
		String start = cols[1];
		String end = cols[2];

		return new Region(Long.parseLong(start), Long.parseLong(end), new Chromosome(chr));
	} 
}
