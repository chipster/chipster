package fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;

import javax.swing.SwingUtilities;

import org.broad.tribble.readers.TabixReader;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Exon;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Gene;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneSet;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.IndexKey;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.BinarySearchIndex;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.GtfLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.Index;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.RandomAccessLineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;

public class GtfToFeatureConversion extends DataThread {
	
	public static int MAX_INTRON_LENGTH = 500*1000; //O,5M should be enough for the longest human introns http://www.bioinfo.de/isb/2004040032/	
	
	private Index index;
	private GtfLineParser parser;
	private boolean isTabix;
	private TabixDataSource tabixDataSource;
	private RandomAccessLineDataSource gtfDataSource;

	public GtfToFeatureConversion(DataUrl gtfTabixUrl, DataUrl gtfIndexUrl, final GBrowser browser) {
	    
		super(browser, null);

		this.isTabix = gtfIndexUrl != null;
		this.parser = new GtfLineParser();
		try {
			
			if (isTabix) {
				tabixDataSource = new TabixDataSource(gtfTabixUrl, gtfIndexUrl);
				super.setDataSource(tabixDataSource);
			} else {
				gtfDataSource = new RandomAccessLineDataSource(gtfTabixUrl);
				this.index = new BinarySearchIndex(gtfDataSource, parser);
				super.setDataSource(gtfDataSource);

				//		try {
				//			this.index = new InMemoryIndex(file, parser);
				//			
				//		} catch (IOException e) {
				//			e.printStackTrace();
				//		}
			}

		} catch (final UnsortedDataException e) {
				SwingUtilities.invokeLater(new Runnable() {

					@Override
					public void run() {				
						browser.showDialog("Unsorted data", e.getMessage(), null, true, false, true, true);
					}
				});
				index = null;

		} catch (IOException e) {
			e.printStackTrace();
		} catch (GBrowserException e) {
			e.printStackTrace();
		} catch (URISyntaxException e) {
			e.printStackTrace();
		}		
	}

	@Override
	protected void processDataRequest(DataRequest request) {				
		
		if (!isTabix && index == null) {
			return;
		}
		

		if (request instanceof GeneRequest) {

			GeneRequest geneRequest = (GeneRequest)request;
			List<Feature> resultList;
			try {
				resultList = processGeneSearch(geneRequest);
				createDataResult(new GeneResult(geneRequest.getStatus(), resultList, geneRequest.getSearchString()));
				
			} catch (IOException e) {
				e.printStackTrace();
			}			

		} else { 

			long start = request.start.bp;
			long end = request.end.bp;

			//Extend area to be able to draw introns at screen edge, but don't below 1
			//TODO Be more clever to avoid getting so much useless data
			start = Math.max((long)start - MAX_INTRON_LENGTH, 1);
			end = end + MAX_INTRON_LENGTH;

			Region requestRegion = new Region(start, end, request.start.chr);

			final long CHUNK_SIZE = 1*1000*1000;

			for (long chunkStart = requestRegion.start.bp; chunkStart <= requestRegion.end.bp; chunkStart += CHUNK_SIZE) {
				Region chunkRegion = new Region(chunkStart, Math.min(chunkStart + CHUNK_SIZE, requestRegion.end.bp), requestRegion.start.chr);

				processDataRequestChunk(request, chunkRegion);
			}
		}
	}
	
	protected void processDataRequestChunk(DataRequest request, Region chunkRegion) {
		
		List<Feature> resultList = new LinkedList<Feature>();
		List<Exon> exons = fetchExons(request, chunkRegion);						
		
		for (Exon exon : exons) {
			
			LinkedHashMap<DataType, Object> valueMap = new LinkedHashMap<DataType, Object>();

			valueMap.put(DataType.VALUE, exon);

			Feature feature = new Feature(exon.getRegion(), valueMap);

			resultList.add(feature);
		}
		
		super.createDataResult(new DataResult(request.getStatus(), resultList));
	}

	private LinkedList<Exon> fetchExons(DataRequest request, Region chunkRegion) {
		
		LinkedList<Exon> exons = new LinkedList<Exon>();
		
		//		long t = System.currentTimeMillis();
		TreeMap<IndexKey, String> lines = null;

		try {
			lines = getLines(request, chunkRegion);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (GBrowserException e) {
			e.printStackTrace();
		}								

		//IndexKeys are not needed, because gtf contains unique identifiers for lines
 		for (String line : lines.values()) {

			parser.setLine(line);					

			Region region = parser.getRegion();
			String feature = parser.getFeature();
			String geneId = parser.getGeneId();
			String transcId = parser.getTranscriptId();

			String exonString = parser.getAttribute("exon_number");
			int exonNumber = -1; 
			if (exonString != null) {
				exonNumber = new Integer(exonString);
			}
			String geneName = parser.getAttribute("gene_name");
			String transcName = parser.getAttribute("transcript_name");
			String biotype = null;

			//Standard gtf data (for example Ensembl)
			if ("exon".equals(feature) || "CDS".equals(feature)) {

				Exon exon = new Exon(region, feature, exonNumber, geneId, transcId, geneName, transcName, biotype);
				exons.add(exon);

				//Custom almost-gtf data
			} else 	if (feature.startsWith("GenBank")) {

				if (geneId == null || transcId == null) {
					continue;
				}

				if ("GenBank gene".equals(feature)) {
					feature = "exon";
				} else if ("GenBank CDS".equals(feature)) {
					feature = "CDS";
				} else {
					geneId = feature + geneId;
					transcId = feature + transcId;

					if (geneName != null) {
						geneName = feature + " " + geneName;
					}

					if (transcName != null) {
						transcName = feature + " " + transcName;
					}

					feature = "exon";
				}

				exonNumber = 1;

				Exon exon = new Exon(region, feature, exonNumber, geneId, transcId, geneName, transcName, biotype);				
				exons.add(exon);
			}
		}
		return exons;
	}

	private TreeMap<IndexKey, String> getLines(DataRequest request,
			Region chunkRegion)
			throws IOException, GBrowserException {
		
		TreeMap<IndexKey, String> lines;
		
		if (isTabix) {
			
			lines = new TreeMap<IndexKey, String>();
			
			TabixReader.Iterator iter = tabixDataSource.getTabixIterator(request);

			String line;
			
			if (iter != null) { //null if there isn't such chromosome in annotations
				
				while ((line = iter.next()) != null) {
					
					parser.setLine(line);
					BpCoord start = parser.getRegion().start;

					
					//create unique indexKeys, otherwise other rows with identical start position are lost
					int i = 0;
					IndexKey key;					
					do {
						key = new IndexKey(start, i);
						i++;
					} while (lines.containsKey(key));
					
					lines.put(key, line);
				}
			}

		} else {		

				lines = index.getFileLines(new DataRequest(chunkRegion, request.getRequestedContents(), request.getStatus()));
		}
		return lines;
	}
	
	private List<Feature> processGeneSearch(GeneRequest request) throws IOException {

		String searchString = request.getSearchString().toLowerCase();
		Chromosome chr = request.start.chr;

		Region region = new Region(1l, Long.MAX_VALUE, chr);

		List<Exon> exons = fetchExons(request, region);
		
		GeneSet genes = new GeneSet();				
		genes.add(exons.iterator(), region);

		List<Feature> resultList = new LinkedList<Feature>();

		for (Gene gene : genes.values()) {

			if (gene.getName() != null && gene.getName().toLowerCase().equals(searchString)) {

				LinkedHashMap<DataType, Object> values = new LinkedHashMap<DataType, Object>();

				values.put(DataType.VALUE, gene);
				resultList.add(new Feature(gene.getRegion(), values));
			}
		}

		return resultList;
	}
}
