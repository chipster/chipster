package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

/**
 * In-memory index, which keeps the whole file in RAM. This is practical for small files
 * below 10 MB, where file reading takes less than 1 second. Files don't need to be sorted. 
 * Memory usage is about 400% in comparison to original file.
 * 
 * @author klemela
 */
public class InMemoryIndex extends Index {

	private LineDataSource file;
	private Parser parser;
	private TreeMap<Region, String> lineMap;

	public InMemoryIndex(DataSource file, Parser parser) throws IOException {
		this.file = (LineDataSource) file;
		this.parser = parser;
		
		readFile();		
	}

	private void readFile() throws IOException {
		
		lineMap = new TreeMap<Region, String>();
		String line;
		while ((line = file.readLine()) != null) {
			
			if (parser.setLine(line)) {
				Region region = parser.getRegion();
				lineMap.put(region, line);
			}
		}
	}
	
	public List<String> getFileLines() {
		
		LinkedList<String> lines = new LinkedList<String>();				
		
		for (Entry<Region, String> entry : lineMap.entrySet()) {

			lines.add(entry.getValue());
		}
		
		return lines;
	}

	public List<String> getFileLines(Region request) {
		
		LinkedList<String> lines = new LinkedList<String>();
		
		Region startRegion = new Region(request.start, request.start);
		Region endRegion = new Region(request.end, request.end);
		
		for (Entry<Region, String> entry : lineMap.subMap(startRegion, endRegion).entrySet()) {

			lines.add(entry.getValue());
		}
		
		return lines;
	}
}
