package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;


public class GeneSearchFileFetcherThread extends Thread {


	private BlockingQueue<BpCoordFileRequest> fileRequestQueue;
	private ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue;

	private GeneSearchHandler areaRequestThread;

	private LineDataSource dataSource;

	private HashMap<String, Chromosome> geneNameMap;

	private  boolean poison = false;

	public GeneSearchFileFetcherThread(BlockingQueue<BpCoordFileRequest> fileRequestQueue, 
			ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue, GeneSearchHandler geneSearchHandler,
			LineDataSource data) {

		this.fileRequestQueue = fileRequestQueue;
		this.fileResultQueue = fileResultQueue;
		this.areaRequestThread = geneSearchHandler;
		this.dataSource = data;

		this.setDaemon(true);
	}

	public void run() {

		while (!poison) {
			try {
				
				for (BpCoordFileRequest fileRequest : fileRequestQueue) {
					if (fileRequest.getStatus().poison) {
						poison = true;
						return;
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

	private void readFile() throws IOException {

		String line;

		while ((line = dataSource.readLine()) != null) {

			String[] cols = line.split("\t");

			if (cols.length == 2) {
				String chr = cols[0];
				String geneName = cols[1];

				geneNameMap.put(geneName.toLowerCase(), new Chromosome(chr));
			}
		}
	}

	private void processFileRequest(BpCoordFileRequest fileRequest) throws IOException {

		if (fileRequest.getStatus().poison) {
			poison = true;
			return;
		}

		if (geneNameMap == null) {
			geneNameMap = new HashMap<String, Chromosome>();

			readFile();
		}

		List<RegionContent> resultList = processGeneSearch((GeneRequest)fileRequest.areaRequest, fileRequest);

		ParsedFileResult result = new ParsedFileResult(resultList, fileRequest, fileRequest.areaRequest, fileRequest.getStatus());

		fileResultQueue.add(result);		
		areaRequestThread.notifyAreaRequestHandler();
	}

	private List<RegionContent> processGeneSearch(GeneRequest areaRequest,
			BpCoordFileRequest fileRequest) throws IOException {

		String searchString = areaRequest.getSearchString();

		Chromosome chr = geneNameMap.get(searchString.toLowerCase());
		List<RegionContent> resultList = new LinkedList<RegionContent>();

		if (chr != null) {

			resultList.add(new RegionContent(new Region(null, null, chr), null));
		}

		return resultList;
	}


	public String toString() {
		return this.getClass().getName() + " - " + dataSource;
	}
}
