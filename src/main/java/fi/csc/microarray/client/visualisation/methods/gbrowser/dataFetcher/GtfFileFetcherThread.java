package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
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
public class GtfFileFetcherThread extends Thread {
	
	private final int COVERAGE_FREQUENCY = 100000;


	private BlockingQueue<BpCoordFileRequest> fileRequestQueue;
	private ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue;

	private GtfHandlerThread areaRequestThread;

	private LineDataSource dataSource;

	private GeneSet genes;
	private TreeMap<Region, Float> coverage = new TreeMap<Region, Float>();
	
	private Chromosome chrInMemory;
	
	private boolean poison = false;


	public GtfFileFetcherThread(BlockingQueue<BpCoordFileRequest> fileRequestQueue, 
			ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue, GtfHandlerThread areaRequestThread,
			LineDataSource dataSource) {

		this.fileRequestQueue = fileRequestQueue;
		this.fileResultQueue = fileResultQueue;
		this.areaRequestThread = areaRequestThread;
		this.dataSource = dataSource;

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
	
	
	private long startTime = 0;
	private long lastTime;
	
	private void stopwatch(String lastOperation) {

//		if (startTime == 0) {
//			startTime = System.currentTimeMillis();
//			lastTime = startTime;
//		} else if ("END".equals(lastOperation)) {
//			System.out.println("TOTAL: \t" + (System.currentTimeMillis() - startTime));
//		} else {
//			System.out.println(lastOperation + ": \t" + (System.currentTimeMillis() - lastTime));
//			lastTime = System.currentTimeMillis();
//		}

	}

	private void readFile() throws IOException {
				
		stopwatch(null);
		
		//Object initialising before loop to avoid unnecessary work for the garbage collecting
		String line;
		
		String[] cols;

		String chr;
		String biotype;
		String feature;
		String exonStart;
		String exonEnd;
		String strand;
		
		String[] ids;

		String geneId;
		String transcId;
		String exonIndex;
		String geneName;
		String transcName;
		
		long skippedRows = 0;
		
		String chrInMemoryString = chrInMemory.toNormalisedString() + "\t";
	
		while ((line = dataSource.readLine()) != null) {
			
			//Check correct chromosome already here to skip uninteresting rows because splitting is slow
			if (!line.startsWith(chrInMemoryString)) {
				continue;
			}
			
			cols = line.split("\t");

			chr = cols[0];
			biotype = cols[1];
			feature = cols[2];
			exonStart = cols[3];
			exonEnd = cols[4];
			strand = cols[6];

			if (biotype.equals("pseudogene")) {
				continue;
			}

			ids = parseIds(cols[8]);
			
			if (ids.length != 5) {
				//Unknown format or missing information
				skippedRows++;
				continue;
			} 
			
			geneId = ids[0];
			transcId = ids[1];
			exonIndex = ids[2];
			geneName = ids[3];
			transcName = ids[4];
			
			Region region = new Region(Long.parseLong(exonStart), Long.parseLong(exonEnd), 
					new Chromosome(chr), getStrand(strand));
		
			Exon exon = new Exon(region, feature, 0); //Integer.parseInt(exonIndex));

			genes.addExon(exon, geneId, transcId, geneName, transcName, biotype);			
		}
		
		if (skippedRows != 0) {
			//TODO replace with proper logging
			System.out.println("Unable to parse " + skippedRows + " rows from file " + dataSource);
		}
		
		stopwatch("read and save");
						
		genes.prepareForReading();
		
		stopwatch("sort");
		
		
		for (Gene gene : genes.values()) {
			for (Transcript transcript : gene.getTranscripts()) {
				for (Exon exon : transcript.getExons()) {
					Region region = exon.getRegion();
					
					long mid = region.getMid();
					Long bin = (mid / COVERAGE_FREQUENCY) * COVERAGE_FREQUENCY; //Rounding to accuracy of COVERAGE_FREQUENCY bp's					
					
					Region binRegion = new Region (bin, bin + COVERAGE_FREQUENCY, 
							exon.getRegion().start.chr.clone(), exon.getRegion().getStrand());
										
					float value = exon.getRegion().getLength();
					
					if (!coverage.containsKey(binRegion)) {
						coverage.put(binRegion, value);
					} else {
						coverage.put(binRegion, coverage.get(binRegion) + value);
					}
					
				}
			}
		}
		
		stopwatch("coverage");
		stopwatch("END");
	} 

	private static Strand getStrand(String strand) {

		if ("+".equals(strand)) {
			return Strand.FORWARD;
		} else if ("-".equals(strand)) {
			return Strand.REVERSED;
		} else if (".".equals(strand)) {
			return Strand.NONE;
		}
		return Strand.UNRECOGNIZED;
	}

	public static String[] parseIds(String ids) {
		String[] split = ids.split(";");

		for (int i = 0; i < split.length; i++) {
			split[i] = split[i].substring(split[i].indexOf("\"") + 1, split[i]
					.lastIndexOf("\""));
		}

		return split;
	}

	private void processFileRequest(BpCoordFileRequest fileRequest) throws IOException {
		
		if (fileRequest.getStatus().poison) {
			poison = true;
			return;
		}

		Chromosome requestChr = fileRequest.areaRequest.start.chr;
		
		if (genes == null || !requestChr.equals(chrInMemory)) {
			genes = new GeneSet();

			chrInMemory = requestChr;
			readFile();
		}

		AreaRequest request = fileRequest.areaRequest;
		List<RegionContent> resultList = new ArrayList<RegionContent>();

		if (!request.status.concise) {

			Collection<Gene> filtered = genes.getGenes(new Region(request.start.bp - 5000, request.end.bp + 5000, new Chromosome("1")));


			for (Gene gene : filtered) {

				LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();

				values.put(ColumnType.VALUE, gene);

				resultList.add(new RegionContent(gene.getRegion(), values));
			}


		} else {
			
			SortedMap<Region, Float> filtered  = coverage.subMap(
					new Region(request.start, request.start), new Region(request.end, request.end));
			

			for (Entry<Region, Float> entry : filtered.entrySet()) {
				
				LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();

				values.put(ColumnType.VALUE, (entry.getValue() / 10000f));
				values.put(ColumnType.STRAND, entry.getKey().getStrand());
				
				resultList.add(new RegionContent(entry.getKey(), values));
			}
		}
		
		ParsedFileResult result = new ParsedFileResult(resultList, fileRequest, request, request.status);

		fileResultQueue.add(result);		
		areaRequestThread.notifyAreaRequestHandler();
	}
	
	public String toString() {
		return this.getClass().getName() + " - " + dataSource;
	}
}
