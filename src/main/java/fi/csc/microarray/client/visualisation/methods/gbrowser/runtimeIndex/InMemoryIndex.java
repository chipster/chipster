package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.IndexKey;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

/**
 * In-memory index, which keeps the whole file in RAM. This is practical for small files
 * below 10 MB, where file reading takes less than 1 second. Files don't need to be sorted. 
 * Memory usage is about 300% in comparison to original file size.
 * 
 * @author klemela
 */
public class InMemoryIndex extends Index {

	private LineDataSource file;
	private LineParser parser;
	
	//TreeMap for storing all lines of the file sorted according to start positions 
	private TreeMap<IndexKey, String> lineMap;

	public InMemoryIndex(DataSource file, LineParser parser) throws IOException {
		this.file = (LineDataSource) file;
		this.parser = parser;
		
		readFile();		
	}

	/**
	 * Read the whole file from disk to RAM
	 * 
	 * @throws IOException
	 */
	private void readFile() throws IOException {
		
		lineMap = new TreeMap<IndexKey, String>();
		String line;
		
		long lineNumber = 0;
		while ((line = file.readLine()) != null) {
			
			if (parser.setLine(line)) {
				IndexKey key = new IndexKey(parser.getRegion().start, lineNumber, true);
				lineMap.put(key, line);
				
				//Count only content lines, because this is a convention in Chipster selection system
				lineNumber++;
			}
			
		}
	}
	
	/**
	 * Get all lines from the file. Use overloaded version {@link #getFileLines(Region)} to get only lines inside specific 
	 * region.
	 * 
	 * @return
	 */
	public List<String> getFileLines() {
		
		LinkedList<String> lines = new LinkedList<String>();				
		
		for (Entry<IndexKey, String> entry : lineMap.entrySet()) {

			lines.add(entry.getValue());
		}
		
		return lines;
	}

	public TreeMap<IndexKey, String> getFileLines(Region request) {
		
		TreeMap<IndexKey, String> lines = new TreeMap<IndexKey, String>();
		
		IndexKey startKey = new IndexKey(request.start, 0l);
		//zero second parameter makes this enKey smaller than lines with equal start position 
		//and thus excludes the lines with equal start position
		IndexKey endKey = new IndexKey(request.end, 0);

		//This shouldn't happen, but lets live on with it
		if (startKey.compareTo(endKey) > 0) {
			System.err.println("Error in class \"" + this.getClass().getSimpleName() + " method getFileLines(): startkey > endKey");
			IndexKey tmp = endKey;
			endKey = startKey;
			startKey = tmp;
		}
		
		for (Entry<IndexKey, String> entry : lineMap.subMap(startKey, endKey).entrySet()) {

			lines.put(entry.getKey(), entry.getValue());
		}
		
		return lines;
	}
	
	@Override
	public Iterator<String> getFileLineIterator() throws IOException, GBrowserException {
		
		return getFileLines().iterator();	
	}
}
