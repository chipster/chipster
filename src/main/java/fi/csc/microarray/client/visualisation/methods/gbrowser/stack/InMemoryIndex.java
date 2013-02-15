package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

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

	public List<String> getFileLines(Region requestRegion) {
		
		LinkedList<String> lines = new LinkedList<String>();				
		
		for (Entry<Region, String> entry : lineMap.entrySet()) {
			
			if (entry.getKey().intersects(requestRegion)) {
				lines.add(entry.getValue());
			}
		}
		
		return lines;
	}
}
