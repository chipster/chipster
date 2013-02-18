package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

public class BinarySearchIndex extends Index {

	private RandomAccessLineDataSource file;
	private Parser parser;
	private TreeMap<BpCoord, Long> index = new TreeMap<BpCoord, Long>();
	
	private static final int INDEX_INTERVAL = 128*1024;

	public BinarySearchIndex(DataSource file, Parser parser) throws IOException, GBrowserException {
		this.file = (RandomAccessLineDataSource) file;
		this.parser = parser;
		
		readEnds();
	}

	private void readEnds() throws IOException, GBrowserException {
		
		//Add first row to index
		file.setLineReaderPosition(0);
		String firstLine = file.getNextLine();		
		parser.setLine(firstLine);
		Region region = parser.getRegion();
		index.put(region.start, 0l);
		
		//Add last row to index
		String lastLine = file.getLastLine();		
		parser.setLine(lastLine);
		region = parser.getRegion();
		index.put(region.start, file.length() - lastLine.length());
	}

	public List<String> getFileLines(Region request) throws IOException, GBrowserException {
		
		//Search start position
		long floorFilePosition = binarySearch(request.start);
		
		file.setLineReaderPosition(floorFilePosition);
		
		LinkedList<String> lines = new LinkedList<String>();
		
		String line = null;
		
		while ((line = file.getNextLine()) != null) {
			
			if ("".equals(line)) {
				//First byte was new line character
				continue;
			}
			
			parser.setLine(line);
			Region region = parser.getRegion();
			
			if (request.contains(region.start)) {
				lines.add(line);
			}
			
			if (request.end.compareTo(region.start) < 0) {
				break;
			}
		}
		
		return lines;
	}
	
	private long binarySearch(BpCoord position) throws IOException, GBrowserException {
		
			Entry<BpCoord, Long> floorEntry = index.floorEntry(position);
			Entry<BpCoord, Long> ceilingEntry = index.ceilingEntry(position);
			
			if (floorEntry == null) {
				//Request start is less than smallest index entry
				floorEntry = ceilingEntry;			
			}
				
			long floorFilePosition = floorEntry.getValue();
			long ceilingFilePosition = ceilingEntry.getValue();
			
			if (ceilingFilePosition - floorFilePosition > INDEX_INTERVAL) {
			
				splitIndexRegion(floorFilePosition, ceilingFilePosition);
				
				return binarySearch(position);
			}
			
			return floorEntry.getValue();
	}

	private void splitIndexRegion(long floorFilePosition,
			long ceilingFilePosition) throws IOException, GBrowserException {
		
		long centerFilePosition = (floorFilePosition + ceilingFilePosition) / 2;
		
		file.setLineReaderPosition(centerFilePosition);
		//Skip the first line, because we don't know if it is complete or partial 		
		long partialLineLength = file.getNextLine().length() + 1; //plus one because of new line character
		String firstLine = file.getNextLine();		
		
		parser.setLine(firstLine);
		Region region = parser.getRegion();
		//File position should be the last new line character before the parsed line.
		//File position cannot be set to actual line start, because then the first line is lost,
		//because we don't know if it is partial or not.
		index.put(region.start, centerFilePosition + partialLineLength + firstLine.length());
	}
}
