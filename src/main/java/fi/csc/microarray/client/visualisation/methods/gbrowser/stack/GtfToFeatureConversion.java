package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.IOException;
import java.util.List;
import java.util.Queue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Exon;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Gene;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneSet;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class GtfToFeatureConversion extends AreaRequestHandler {

	private InMemoryIndex index;

	private StackGtfParser parser;

	public GtfToFeatureConversion(DataSource file, Queue<AreaRequest> areaRequestQueue,
	        AreaResultListener areaResultListener) {
	    
		super(areaRequestQueue, areaResultListener);

		this.parser = new StackGtfParser();
		
		try {
			this.index = new InMemoryIndex(file, parser);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Override
	protected void processAreaRequest(AreaRequest request) {
						
		super.processAreaRequest(request);	
		
		if (request.getStatus().poison) {
			return;
		}
		
		//limit to integer range
		long start = request.start.bp;
		long end = request.end.bp;
		
		//Extend area to be able to draw introns at screen edge, but don't go over MAX_VALUE, or below 1
		//TODO Be more clever to avoid getting so much useless data
		int EXTRA = 500000; //O,5M should be enough for the longest human introns http://www.bioinfo.de/isb/2004040032/
		
		start = Math.max((long)start - EXTRA, 1);
		end = end + EXTRA;
		
		Region requestRegion = new Region(start, end, request.start.chr);
		List<String> lines = index.getFileLines(new AreaRequest(requestRegion, request.getRequestedContents(), request.getStatus()));
		
		GeneSet geneSet = new GeneSet();
		
		for (String line : lines) {
			
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
		
		for (Gene gene : geneSet.values()) {
						
			super.createAreaResult(new AreaResult(request.getStatus(), gene.getRegion(), gene));
		}
	}
}
