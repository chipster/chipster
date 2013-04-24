package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.IndexKey;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;

/**
 * This is a general conversion class for tracks that need only the region information. All file parsing
 * should happen here, and the data should be converted so that the track drawing is as fast as possible.
 * The results include the original line, but for performance reasons it should be used only for debug 
 * etc. purposes.
 *  
 * Conversion class that extracts region information from the lines using the supplied parser.
 * Results are in RegionContent objects. These objects contain:
 * <ul>
 * <li>the region
 * <li>an unique line identifier stored with key ColumnType.ID
 * <li>the original line stored with key ColumnType.VALUE
 * </ul> 
 * 
 * @author klemela
 *
 */
public class LineToRegionConversion extends SingleThreadAreaRequestHandler {

	private Index index;

	private LineParser parser;

	public LineToRegionConversion(DataSource file, LineParser parser, Queue<AreaRequest> areaRequestQueue,
	        AreaResultListener areaResultListener) {
	    
		super(areaRequestQueue, areaResultListener);

		this.parser = parser;
		
//		try {
//			this.index = new InMemoryIndex(file, parser);
//			
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
		try {
			this.index = new BinarySearchIndex(file, parser);
		
		} catch (UnsortedDataException e) {
			e.printStackTrace();
			index = null;
		} catch (IOException e) {
			e.printStackTrace();
		} catch (GBrowserException e) {
			e.printStackTrace();
		}
	}

	public LineToRegionConversion(URL url, LineParser parser) throws FileNotFoundException, URISyntaxException {
		this(new RandomAccessLineDataSource(url), parser, null, null);
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
		int EXTRA = 1000000; 
		
		start = Math.max((long)start - EXTRA, 1);
		end = end + EXTRA;
		
		Region requestRegion = new Region(start, end, request.start.chr);
				
		TreeMap<IndexKey, String> lines = null;
		try {		
					
			lines = index.getFileLines(new AreaRequest(requestRegion, request.getRequestedContents(), request.getStatus()));
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (GBrowserException e) {
			e.printStackTrace();
		}	
				
		List<RegionContent> list = new LinkedList<RegionContent>();
		
		for (Entry<IndexKey, String> entry : lines.entrySet()) {
			
			String line = entry.getValue();			
			parser.setLine(line);								
			Region region = parser.getRegion();
			IndexKey id = entry.getKey();
			
			LinkedHashMap<ColumnType, Object> valueMap = new LinkedHashMap<ColumnType, Object>();			
			
			valueMap.put(ColumnType.ID, id);
			valueMap.put(ColumnType.VALUE, line);			
			RegionContent regionContent = new RegionContent(region, valueMap);
			
			list.add(regionContent);
		}						
			
		super.createAreaResult(new AreaResult(request.getStatus(), list));
	}
}
