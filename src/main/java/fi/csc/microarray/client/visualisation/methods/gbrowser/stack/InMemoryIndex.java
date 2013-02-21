package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

/**
 * In-memory index, which keeps the whole file in RAM. This is practical for small files
 * below 10 MB, where file reading takes less than 1 second. Files don't need to be sorted. 
 * Memory usage is about 300% in comparison to original file.
 * 
 * Example output:
 * 
 * Memory usage: 0 MB
 * Init memIndex: 24705 ms
 * Memory usage: 1528 MB
 * Init  searchIndex: 10 ms
 * Memory usage: 1528 MB
 * Troughput test
 * memIndex:     	Bytes: 20 MB 	Time: 54 ms 	Bandwidth: 377.79705612747756MB/s
 * searchIndex:  	Bytes: 20 MB 	Time: 3575 ms 	Bandwidth: 5.706584903743717MB/s
 * Random access test
 * memIndex:     Seek 1000, first: 18 ms
 * Memory usage: 1573 MB
 * memIndex:    Seek 1000, second: 11 ms
 * Memory usage: 1576 MB
 * searchIndex:  Seek 1000, first: 13673 ms
 * Memory usage: 1587 MB
 * searchIndex: Seek 1000, second: 13862 ms
 * Memory usage: 1540 MB
 * Agreement test
 * Agreement test done
 * 
 * @author klemela
 */
public class InMemoryIndex extends Index {
	
	private class IndexKey  implements Comparable<IndexKey> {
		
		public IndexKey(BpCoord start, long lineNumber) {
			this.start = start;
			this.lineNumber = lineNumber;
		}
		
		private BpCoord start;
		private long lineNumber;

		@Override
		public boolean equals(Object o) {
			if (o instanceof IndexKey) {
				return (this.compareTo((IndexKey) o) == 0);
			}
			return false;
		}

		@Override
		public int hashCode() {
			return start.hashCode();
		}

		@Override
		public int compareTo(IndexKey o) {

			int startComparison = start.compareTo(o.start);

			if (startComparison != 0) {
				return startComparison;

			} else {
				return ((Long)lineNumber).compareTo(o.lineNumber);
			}
		}
	}

	private LineDataSource file;
	private Parser parser;
	private TreeMap<IndexKey, String> lineMap;

	public InMemoryIndex(DataSource file, Parser parser) throws IOException {
		this.file = (LineDataSource) file;
		this.parser = parser;
		
		readFile();		
	}

	private void readFile() throws IOException {
		
		lineMap = new TreeMap<IndexKey, String>();
		String line;
		
		long lineNumber = 0;
		while ((line = file.readLine()) != null) {
			
			if (parser.setLine(line)) {
				IndexKey key = new IndexKey(parser.getRegion().start, lineNumber);
				lineMap.put(key, line);				
			}
			
			lineNumber++;
		}
	}
	
	public List<String> getFileLines() {
		
		LinkedList<String> lines = new LinkedList<String>();				
		
		for (Entry<IndexKey, String> entry : lineMap.entrySet()) {

			lines.add(entry.getValue());
		}
		
		return lines;
	}

	public List<String> getFileLines(Region request) {
		
		LinkedList<String> lines = new LinkedList<String>();
		
		IndexKey startKey = new IndexKey(request.start, 0l);
		IndexKey endKey = new IndexKey(request.end, Long.MAX_VALUE);
		
		for (Entry<IndexKey, String> entry : lineMap.subMap(startKey, endKey).entrySet()) {

			lines.add(entry.getValue());
		}
		
		return lines;
	}
}
