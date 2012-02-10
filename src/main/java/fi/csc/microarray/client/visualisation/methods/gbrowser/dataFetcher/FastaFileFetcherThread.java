package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.FastaDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * 
 * The data retrieval layer thread for gtf files from the ftp.ensembl.org containing genes and transcripts. 
 * Keeps annotations of only one chromosome in memory to keep memory consumption in about 200 MB (for 1. human 
 * chromosome). Annotation file is about 400 MB and has to be read through in every chromosome change, so
 * performance and memory consumption of readFile() method is quite important. This allows showing annotations
 * without too much waiting, currently in about 4 seconds. Receives file requests and sends file results.
 * 
 * Data is stored in hierarchical structure of Gene, Transcript and Exon objects. After file reading a 
 * coverage is calculated representing number of exons over the chromosome.
 * 
 * @author Petri Klemelä
 *
 */
public class FastaFileFetcherThread extends Thread {

	private BlockingQueue<BpCoordFileRequest> fileRequestQueue;
	private ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue;

	private FastaHandlerThread areaRequestThread;
	
	private Map<Chromosome, Fasta> fastas;
	
	private boolean poison = false;
	
	
	private class Fasta {
		
		public Fasta(ChunkDataSource dataSource) {
			this.dataSource = dataSource;
		}
		
		public void init() throws IOException {
			
			byte[] headBytes = new byte[10*1000]; //10k max total length of header and first data line
			
			dataSource.read(0, headBytes);
			
			String headString = new String(headBytes);
			
			int firstNewLine = headString.indexOf("\n");
			int secondNewLine = headString.indexOf("\n", firstNewLine + 1);
			
			headerLength = firstNewLine + 1; //plus one because of new line character
			fileRowLength = secondNewLine - firstNewLine; //including new line character
		}
		
		public long bpToFile(long bp) {
			long bpRowLength = fileRowLength - 1;//Minus one because file rowLength contains new line character
			int rows = (int) (bp / bpRowLength);//Minus one to calculate in 0-based coordinate system
			int column = (int) (bp % bpRowLength); //Minus one to calculate in 0-based coordinate system
			return headerLength + rows * fileRowLength + column;
		}
		
		public String read(long bpStart, long bpEnd) throws IOException {
			
			//Convert to 0-based coordinate system
			bpStart--;
			bpEnd--;
			
			if (headerLength == -1) {
				init();
			}
			
			long startPosition = bpToFile(bpStart);
			long endPosition = bpToFile(bpEnd);
			
			byte[] bytes = new byte[(int) (endPosition - startPosition) + 1];
			
			dataSource.read(startPosition, bytes);
			
			String rowSepareated = new String(bytes);
			
			String requestedData = rowSepareated.replace("\n", "");
			
			return requestedData;
		}
		
		private ChunkDataSource dataSource;
		private long headerLength = -1;
		private long fileRowLength = -1;
	}

	public FastaFileFetcherThread(BlockingQueue<BpCoordFileRequest> fileRequestQueue, 
			ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue, FastaHandlerThread areaRequestThread,
			FastaDataSource dataSource) {

		this.fileRequestQueue = fileRequestQueue;
		this.fileResultQueue = fileResultQueue;
		this.areaRequestThread = areaRequestThread;
		
		this.fastas = new TreeMap<Chromosome, Fasta>();
		
		for (Entry<Chromosome, ChunkDataSource> entry : dataSource.entrySet()) {
			fastas.put(entry.getKey(), new Fasta(entry.getValue()));
		}

		this.setDaemon(true);
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
	}

	private void processFileRequest(BpCoordFileRequest fileRequest) throws IOException {

		if (fileRequest.getStatus().poison) {
			poison = true;
			return;
		}

		AreaRequest request = fileRequest.areaRequest;

		List<RegionContent> resultList = new ArrayList<RegionContent>();

		Chromosome chr = request.start.chr;
		
		Fasta fasta = fastas.get(chr);

		String seqence = fasta.read(request.start.bp, request.end.bp);
		
		
		LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();
		
		values.put(ColumnType.SEQUENCE, seqence);

		resultList.add(new RegionContent(new Region(request), values));


		ParsedFileResult result = new ParsedFileResult(resultList, fileRequest, request, request.status);

		fileResultQueue.add(result);		
		areaRequestThread.notifyAreaRequestHandler();
	}
	
	public String toString() {
		return this.getClass().getName() + " - " + fastas;
	}
}
