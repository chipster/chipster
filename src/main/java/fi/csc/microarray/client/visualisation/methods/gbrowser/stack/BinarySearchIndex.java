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
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;

public class BinarySearchIndex extends Index {

	private RandomAccessLineDataSource file;
	private Parser parser;
	private TreeMap<BpCoord, Long> index = new TreeMap<BpCoord, Long>();
	
	private static final int INDEX_INTERVAL = 128*1024;

	public BinarySearchIndex(DataSource file, Parser parser) throws IOException, GBrowserException {
		this.file = (RandomAccessLineDataSource) file;
		this.parser = parser;
		
		checkSorting();
		
		readEnds();
	}

	public void checkSorting() throws IOException, GBrowserException {

		List<String> lines;

		if (getFile().length() < 100*1024) {
			lines = getFileLines();

		} else {

			long length = getFile().length();

			final int blockLineCount = 10;
			final int blockCount = 10;

			final long blockInterval = length / blockCount;

			lines = new LinkedList<String>();

			for (int i = 0; i < blockCount; i++) {

				getFile().setLineReaderPosition(blockInterval * i);
				getFile().getNextLine(); //Maybe partial line

				for (int j = 0; j < blockLineCount; j++) {
					lines.add(getFile().getNextLine());
				}
			}
		}		

		checkSorting(lines);
	}

	private void checkSorting(List<String> lines) throws UnsortedDataException {
		
		String previousLine = null;
		Region previousRegion = null;
		
		for (String line : lines) {
			
			getParser().setLine(line);
			Region region = getParser().getRegion();
			
			if (previousRegion != null) {
				if (previousRegion.start.compareTo(region.start) > 0) {
					throw new UnsortedDataException("File " + getFile() + " isn't sorted correctly. " +
							"Please sort the file first. First line: " + previousLine + " Second line: " + line);
				}
			}
			previousRegion = region;
			previousLine = line;
		}	
	}

	private void readEnds() throws IOException, GBrowserException {
		
		//Add first row to index
		getFile().setLineReaderPosition(0);
		String firstLine = getFile().getNextLine();		
		getParser().setLine(firstLine);
		Region region = getParser().getRegion();
		index.put(region.start, 0l);
		
		//Add last row to index
		String lastLine = getFile().getLastLine();		
		getParser().setLine(lastLine);
		region = getParser().getRegion();
		index.put(region.start, getFile().length() - lastLine.length());
	}
	
	/**
	 * Read the whole file and return a list of lines. Obviously, use this only for small
	 * files.
	 * 
	 * @return
	 * @throws IOException
	 * @throws GBrowserException
	 */
	public List<String> getFileLines() throws IOException, GBrowserException {
				
		getFile().setLineReaderPosition(0);
		
		LinkedList<String> lines = new LinkedList<String>();
		
		String line = null;
		
		while ((line = getFile().getNextLine()) != null) {
						
			getParser().setLine(line);
			Region region = getParser().getRegion();
			
			lines.add(line);
		}
		
		return lines;
	}

	public List<String> getFileLines(Region request) throws IOException, GBrowserException {
		
		//Search start position
		long floorFilePosition = binarySearch(request.start);
		
		getFile().setLineReaderPosition(floorFilePosition);
		
		LinkedList<String> lines = new LinkedList<String>();
		
		String line = null;
		
		while ((line = getFile().getNextLine()) != null) {
			
			if ("".equals(line)) {
				//First byte was new line character
				continue;
			}
			
			getParser().setLine(line);
			Region region = getParser().getRegion();
			
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
			
			if (ceilingEntry == null) {
				//Request start is greater than largest index entry
				ceilingEntry = floorEntry;			
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
		
		getFile().setLineReaderPosition(centerFilePosition);
		//Skip the first line, because we don't know if it is complete or partial 		
		long partialLineLength = getFile().getNextLine().length() + 1; //plus one because of new line character
		String firstLine = getFile().getNextLine();		
		
		getParser().setLine(firstLine);
		Region region = getParser().getRegion();
		//File position should be the last new line character before the parsed line.
		//File position cannot be set to actual line start, because then the first line is lost,
		//because we don't know if it is partial or not.
		index.put(region.start, centerFilePosition + partialLineLength + firstLine.length());
	}

	public RandomAccessLineDataSource getFile() {
		return file;
	}

	public Parser getParser() {
		return parser;
	}
}
