package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;

import javax.swing.SwingUtilities;

import org.broad.tribble.readers.TabixReader;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Exon;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.TabixUtil;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;

public class GtfToFeatureConversion extends SingleThreadAreaRequestHandler {
	
	public static int MAX_INTRON_LENGTH = 500*1000; //O,5M should be enough for the longest human introns http://www.bioinfo.de/isb/2004040032/	
	
	private Index index;
	private GtfLineParser parser;
	private boolean isTabix;
	private TabixDataSource tabixDataSource;
	private RandomAccessLineDataSource gtfDataSource;

	public GtfToFeatureConversion(URL dataUrl, URL indexUrl, final GBrowser browser) {
	    
		super(null, null);

		this.isTabix = indexUrl != null;
		this.parser = new GtfLineParser();
		try {
			
			if (isTabix) {
				tabixDataSource = new TabixDataSource(dataUrl, indexUrl);
			} else {
				gtfDataSource = new RandomAccessLineDataSource(dataUrl);
				this.index = new BinarySearchIndex(gtfDataSource, parser);

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
	protected void processAreaRequest(AreaRequest request) {
						
		super.processAreaRequest(request);	
		
		if (request.getStatus().poison) {
			return;
		}
		
		if (!isTabix && index == null) {
			return;
		}
		
		long start = request.start.bp;
		long end = request.end.bp;
		
		//Extend area to be able to draw introns at screen edge, but don't below 1
		//TODO Be more clever to avoid getting so much useless data
		start = Math.max((long)start - MAX_INTRON_LENGTH, 1);
		end = end + MAX_INTRON_LENGTH;
		
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
			lines = getLines(request, chunkRegion);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (GBrowserException e) {
			e.printStackTrace();
		}
		
		List<RegionContent> list = new LinkedList<RegionContent>();
		
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
			
			//Standard gtf data (for example Ensembl)
			if ("exon".equals(feature) || "CDS".equals(feature)) {

				Exon exon = new Exon(region, feature, exonNumber, geneId, transcId, geneName, transcName, biotype);				
				addExon(list, exon);
				
			//Custom almost-gtf data
			} else 	if (feature.startsWith("GenBank")) {
				
				if (geneId == null || transcId == null) {
					return;
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
				addExon(list, exon);
			}
		}						
		
		super.createAreaResult(new AreaResult(request.getStatus(), list));
	}
	
	private void addExon(List<RegionContent> list, Exon exon) {
		
		LinkedHashMap<ColumnType, Object> valueMap = new LinkedHashMap<ColumnType, Object>();
		
		valueMap.put(ColumnType.VALUE, exon);
		
		RegionContent feature = new RegionContent(exon.getRegion(), valueMap);
		
		list.add(feature);
	}

	private TreeMap<IndexKey, String> getLines(AreaRequest request,
			Region chunkRegion)
			throws IOException, GBrowserException {
		
		TreeMap<IndexKey, String> lines;
		
		if (isTabix) {
			
			lines = new TreeMap<IndexKey, String>();
			
			TabixReader.Iterator iter = TabixUtil.getTabixIterator(tabixDataSource, request);

			String line;
			
			if (iter != null) { //null if there isn't such chromosome in annotations
				
				while ((line = iter.next()) != null) {
					
					parser.setLine(line);
					BpCoord start = parser.getRegion().start;
					//proper IndexKeys are not needed, because gtf contains unique identifiers for lines
					lines.put(new IndexKey(start, 0), line);
				}
			}

		} else {		

				lines = index.getFileLines(new AreaRequest(chunkRegion, request.getRequestedContents(), request.getStatus()));
		}
		return lines;
	}
}
