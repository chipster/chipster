package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.SortedMap;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class BinarySearchIndex extends Index {

	private ChunkDataSource file;
	private Parser parser;
	private TreeMap<BpCoord, Long> index = new TreeMap<BpCoord, Long>();
	
	private static final int CHUNK_LENGTH = 1024;
	private static final int INDEX_INTERVAL = 128;

	public BinarySearchIndex(DataSource file, Parser parser) throws IOException {
		this.file = (ChunkDataSource) file;
		this.parser = parser;
		
		readEnds();
	}

	private void readEnds() throws IOException {
		byte[] chunk = new byte[CHUNK_LENGTH];
		
		//Add first row to index
		file.read(0, chunk);
		String chunkString = new String(chunk);
		String firstLine = chunkString.substring(0, chunkString.indexOf("\n"));
		
		parser.setLine(firstLine);
		Region region = parser.getRegion();
		index.put(region.start, 0l);
		
		//Add last row to index
		long chunkPosition = file.length() - CHUNK_LENGTH;
		file.read(chunkPosition, chunk);
		chunkString = new String(chunk);
		chunkString = chunkString.substring(0, chunkString.lastIndexOf("\n"));
		int indexOfLastLine = chunkString.lastIndexOf("\n") + 1;
		String lastLine = chunkString.substring(indexOfLastLine);
		
		parser.setLine(lastLine);
		region = parser.getRegion();
		index.put(region.start, chunkPosition + indexOfLastLine);
	}

	public List<String> getFileLines(Region request) {
		
		//Update index around request start
		binarySearch(request.start);
		
		Entry<BpCoord, Long> floorEntry = index.floorEntry(request.start);
		
		long floorFilePosition = floorEntry.getValue();
		
		//Get requested region from index
		SortedMap<BpCoord, Long> requestIndex = index.subMap(request.start, request.end);
		
		LinkedList<String> lines = new LinkedList<String>();				
		
		for (Entry<Region, String> entry : lineMap.entrySet()) {
			
			if (entry.getKey().intersects(requestRegion)) {
				lines.add(entry.getValue());
			}
		}
		
		return lines;
	}
	
	private long binarySearch(BpCoord position) throws IOException {
		
			Entry<BpCoord, Long> floorEntry = index.floorEntry(position);
			Entry<BpCoord, Long> ceilingEntry = index.ceilingEntry(position);
			
			long floorFilePosition = floorEntry.getValue();
			long ceilingFilePosition = ceilingEntry.getValue();
			
			if (ceilingFilePosition - floorFilePosition > CHUNK_LENGTH * INDEX_INTERVAL) {
			
				splitIndexRegion(floorFilePosition, ceilingFilePosition);
				
				return binarySearch(position);
			}
			
			return floorEntry.getValue();
	}

	private void splitIndexRegion(long floorFilePosition,
			long ceilingFilePosition) throws IOException {
		
		long centerFilePosition = (floorFilePosition + ceilingFilePosition) / 2;
		
		byte[] chunk = new byte[CHUNK_LENGTH];

		file.read(centerFilePosition, chunk);
		String chunkString = new String(chunk);
		int indexOfFirstLine = chunkString.indexOf("\n") + 1;
		chunkString = chunkString.substring(indexOfFirstLine);
		String firstLine = chunkString.substring(0, chunkString.lastIndexOf("\n"));
		
		parser.setLine(firstLine);
		Region region = parser.getRegion();
		index.put(region.start, centerFilePosition + indexOfFirstLine);
	}
}
