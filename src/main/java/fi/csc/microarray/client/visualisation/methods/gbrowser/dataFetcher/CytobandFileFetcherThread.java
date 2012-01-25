package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.CytobandDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * 
 * The data retrieval layer thread for Cytoband files from the ftp.ensembl.org . Uses the main karyotype file 
 * and region file to interpret region codes to chromosomes. Keeps everything in memory, which should not 
 * be problem with these small files. Receives file requests and sends file results.
 * 
 * @author Petri Klemelä
 *
 */
public class CytobandFileFetcherThread extends Thread {


	private BlockingQueue<BpCoordFileRequest> fileRequestQueue;
	private ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue;
	
	private CytobandHandlerThread areaRequestThread;

	private CytobandDataSource dataSource;
		
	SortedSet<Cytoband> cytobands;

	public CytobandFileFetcherThread(BlockingQueue<BpCoordFileRequest> fileRequestQueue, 
			ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue, CytobandHandlerThread areaRequestThread,
			CytobandDataSource dataSource) {

		this.fileRequestQueue = fileRequestQueue;
		this.fileResultQueue = fileResultQueue;
		this.areaRequestThread = areaRequestThread;
		this.dataSource = dataSource;

		this.setDaemon(true);
	}

	public void run() {

		while (true) {
			try {
				processFileRequest(fileRequestQueue.take());

			} catch (IOException e) {
				e.printStackTrace(); // FIXME fix exception handling
			} catch (InterruptedException e) {
				e.printStackTrace(); // FIXME fix exception handling
			}
		}
	}

	private void readFile() throws IOException {
		
				Map<String, String> chrMap = new HashMap<String, String>();
				
				String line;
				
				LineDataSource regionDataSource = dataSource.getRegionDataSrouce();
				
				while ((line = regionDataSource.readLine()) != null) {

					String[] cols = line.split("\t");

					//See ensembl mysql schema
					String seq_region_id = cols[0];
					String chr = cols[1];
					String coord_system_id = cols[2];
					//String length = cols[3];
					
					
					//seq_region types can be found from following file, but now it's assumed that type 2 means
					//chromosome: ftp://ftp.ensembl.org/pub/release-57/mysql/homo_sapiens_core_57_37b/coord_system.txt.gz
					if (!(coord_system_id.equals("2") || coord_system_id.equals("101"))) {
						continue;
					}
										
					chrMap.put(seq_region_id, chr);
				}
				
				while ((line = dataSource.readLine()) != null) {

					String[] cols = line.split("\t");

					//See ensembl mysql schema
					//String karyotype_id = cols[0];
					String seq_region_id = cols[1];
					String seq_region_start = cols[2];
					String seq_region_end = cols[3];
					String band = cols[4];
					String stain = cols[5];
					
					String chr = chrMap.get(seq_region_id);
					
					if (chr == null) {

						continue;
						
					} else {

						BpCoordRegion region = new BpCoordRegion(
								Long.decode(seq_region_start), Long.decode(seq_region_end), new Chromosome(chr));

						Cytoband cytoband = new Cytoband(region, band, stain);

						cytobands.add(cytoband);
					}
				}
	}

	private void processFileRequest(BpCoordFileRequest fileRequest) throws IOException {
		
		if (cytobands == null) {
			cytobands = new TreeSet<Cytoband>();
			
			readFile();
		}
			
		AreaRequest request = fileRequest.areaRequest;

		SortedSet<Cytoband> filtered = cytobands.subSet(
				new Cytoband(request.start), new Cytoband(request.end));
		
		List<RegionContent> resultList = new ArrayList<RegionContent>();
		
		
		for (Cytoband cband : filtered) {
			
			LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();
			
			values.put(ColumnType.VALUE, cband);
			
			resultList.add(new RegionContent(cband.getRegion(), values));
		}
			
		ParsedFileResult result = new ParsedFileResult(resultList, fileRequest, request, request.status);
				
		fileResultQueue.add(result);		
		areaRequestThread.notifyAreaRequestHandler();
	}
}
