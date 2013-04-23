package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Exon;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Gene;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneSet;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;

public class GtfToFeatureConversion extends SingleThreadAreaRequestHandler {

	private Index index;

	private GtfLineParser parser;

	public GtfToFeatureConversion(DataSource file, final GBrowser browser) {
	    
		super(null, null);

		this.parser = new GtfLineParser();
		
//		try {
//			this.index = new InMemoryIndex(file, parser);
//			
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
		try {
			this.index = new BinarySearchIndex(file, parser);
		
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
		}
	}

	@Override
	protected void processAreaRequest(AreaRequest request) {
						
		super.processAreaRequest(request);	
		
		if (request.getStatus().poison) {
			return;
		}
		
		if (index == null) {
			return;
		}
		
		long start = request.start.bp;
		long end = request.end.bp;
		
		//Extend area to be able to draw introns at screen edge, but don't below 1
		//TODO Be more clever to avoid getting so much useless data
		int EXTRA = 500000; //O,5M should be enough for the longest human introns http://www.bioinfo.de/isb/2004040032/
		
		start = Math.max((long)start - EXTRA, 1);
		end = end + EXTRA;
		
		Region requestRegion = new Region(start, end, request.start.chr);
		
		final long CHUNK_SIZE = 1*1000*1000;
		if (requestRegion.getLength() > CHUNK_SIZE) {
		
			for (long chunkStart = requestRegion.start.bp; chunkStart <= requestRegion.end.bp; chunkStart += CHUNK_SIZE) {
				Region chunkRegion = new Region(chunkStart, Math.min(chunkStart + CHUNK_SIZE, requestRegion.end.bp), requestRegion.start.chr);
				
				processAreaRequestChunk(request, chunkRegion);
			}
		}		
	}
	
	protected void processAreaRequestChunk(AreaRequest request, Region chunkRegion) {
		
//		long t = System.currentTimeMillis();
		
		TreeMap<IndexKey, String> lines = null;
		try {		
			
			lines = index.getFileLines(new AreaRequest(chunkRegion, request.getRequestedContents(), request.getStatus()));
			
//			System.out.println("getFileLines\t " + (System.currentTimeMillis() - t) + " ms, lines:\t " + lines.size());
//			t = System.currentTimeMillis();
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (GBrowserException e) {
			e.printStackTrace();
		}
		
		GeneSet geneSet = new GeneSet();
		
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
			String transcName = parser.getAttribute("gene_name");
			String geneName = parser.getAttribute("transcript_name");
			String biotype = null;
			
			Exon exon = new Exon(region, feature, exonNumber);
			
			geneSet.addExon(exon, geneId, transcId, geneName, transcName, biotype);
		}					
		
		List<RegionContent> list = new LinkedList<RegionContent>();
		
		for (Gene gene : geneSet.values()) {
						
			LinkedHashMap<ColumnType, Object> valueMap = new LinkedHashMap<ColumnType, Object>();			
			valueMap.put(ColumnType.VALUE, gene);			
			RegionContent regionContent = new RegionContent(gene.getRegion(), valueMap);
			
			list.add(regionContent);
		}		
		
		super.createAreaResult(new AreaResult(request.getStatus(), list));
	}
}
