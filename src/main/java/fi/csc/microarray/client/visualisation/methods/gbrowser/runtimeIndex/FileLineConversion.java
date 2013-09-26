package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.IndexKey;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;

/**
 * This class converts bed files to BedLine objects.
 *  
 * Results are in Feature objects. These objects contain:
 * <ul>
 * <li>the region
 * <li>an unique line identifier stored with key DataType.ID
 * <li>the BedLine object stored with key DataType.VALUE
 * </ul> 
 * 
 * @author klemela
 *
 */
public class FileLineConversion extends DataThread {

	//Use InMemoryIndex for files under 1MB
	private static final long IN_MEMORY_INDEX_LIMIT = 1*1000*1000;

	private Index index;

	private LineParser parser;

	public FileLineConversion(DataUrl data, AbstractTsvLineParser parser, GBrowser browser) throws FileNotFoundException, URISyntaxException {
		super(browser, null);
						
		try {					
						 			
			RandomAccessLineDataSource dataSource = new RandomAccessLineDataSource(data);
			this.parser = parser;
			
			if (dataSource.length() < IN_MEMORY_INDEX_LIMIT) {
				//InMemoryIndex requires different kind of DataSource
				LineDataSource lineDataSource = new LineDataSource(data);
				this.index = new InMemoryIndex(lineDataSource, parser);
				super.setDataSource(lineDataSource);
			} else {
				this.index = new BinarySearchIndex(dataSource, parser);
				super.setDataSource(dataSource);
			}
		
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
				
		List<Feature> list = new LinkedList<Feature>();
		
		for (Entry<IndexKey, String> entry : lines.entrySet()) {
			
			String line = entry.getValue();			
			parser.setLine(line);								
			Region region = parser.getRegion();
			FileLine fileLine = parser.getFileLine();
			IndexKey id = entry.getKey();
			
			LinkedHashMap<DataType, Object> valueMap = new LinkedHashMap<DataType, Object>();			
			
			valueMap.put(DataType.ID, id);
			valueMap.put(DataType.VALUE, fileLine);			
			Feature regionContent = new Feature(region, valueMap);
			
			list.add(regionContent);
		}						
			
		super.createDataResult(new DataResult(request.getStatus(), list));
	}

	public LineParser getParser() {
		return parser;
	}

	public Index getIndex() {
		return index;
	}
}
