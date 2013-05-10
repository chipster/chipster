package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
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
 * <li>an unique line identifier stored with key DataType.ID
 * <li>the original line stored with key DataType.VALUE
 * </ul> 
 * 
 * @author klemela
 *
 */
public class LineToRegionConversion extends DataThread {

	private Index index;

	private LineParser parser;

	public LineToRegionConversion(URL url, LineParser parser, GBrowser browser) throws FileNotFoundException, URISyntaxException {
		super(browser);
			    
		this.parser = parser;
		
//		try {
//			this.index = new InMemoryIndex(file, parser);
//			
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
		try {
			this.index = new BinarySearchIndex(new RandomAccessLineDataSource(url), parser);
		
		} catch (UnsortedDataException e) {
			e.printStackTrace();
			index = null;
		} catch (IOException e) {
			e.printStackTrace();
		} catch (GBrowserException e) {
			e.printStackTrace();
		}
	}

	@Override
	protected void processDataRequest(DataRequest request) {					
		
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
					
			lines = index.getFileLines(new DataRequest(requestRegion, request.getRequestedContents(), request.getStatus()));
			
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
			
			LinkedHashMap<DataType, Object> valueMap = new LinkedHashMap<DataType, Object>();			
			
			valueMap.put(DataType.ID, id);
			valueMap.put(DataType.VALUE, line);			
			RegionContent regionContent = new RegionContent(region, valueMap);
			
			list.add(regionContent);
		}						
			
		super.createDataResult(new DataResult(request.getStatus(), list));
	}
}
