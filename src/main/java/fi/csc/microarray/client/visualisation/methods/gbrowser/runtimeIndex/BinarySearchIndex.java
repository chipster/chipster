package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Map.Entry;
import java.util.TreeMap;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.IndexKey;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;

/**
 * BinarySearchIndex locates requested region from the sorted file. The requested region is 
 * located by binary search algorithm. First searches require disk seeks, but most subsequent 
 * searches can be satisfied with RAM index, which is build dynamically during the searches.
 * The actual content is read always from the file, keeping the memory usage minimal.
 * 
 * Comment and header lines are tolerated only in the beginning of the file.
 * 
 * @author klemela
 */
public class BinarySearchIndex extends Index {

	private RandomAccessLineDataSource file;
	private LineParser parser; 
	
	/**
	 * Index maps BpCoordinates (i.e. the value of start column of line) to file position of the
	 * preceding new line character (in bytes).  
	 *  
	 * The file positions in the index don't point to the first character of the line, but to the 
	 * last new line character before it. When that location is read from file, the new line character 
	 * is interpreted to mean that next line is complete.
	 */
	private TreeMap<BpCoord, Long> index = new TreeMap<BpCoord, Long>();
	
	private static final int INDEX_INTERVAL = 128*1024;

	public BinarySearchIndex(DataSource file, LineParser parser) throws IOException, GBrowserException, UnsortedDataException {
		this.file = (RandomAccessLineDataSource) file;
		this.parser = parser;
		
		checkSorting();		
		readEnds();
	}

	/**
	 * Check that file is sorted correctly. Small files are checked completely, but
	 * in bigger files only 10 chunks of 10 lines are picked for this sorting check.
	 *  
	 * @throws UnsortedDataException Is thrown if the sorting isn't correct. 
	 * @throws IOException
	 * @throws GBrowserException
	 */
	public void checkSorting() throws IOException, GBrowserException, UnsortedDataException {

		List<String> lines;

		// the limit should be greater than blockLineCount * blockCount * maxLineLength to avoid
		// problems in big file implementation
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
					String line = getFile().getNextLine();
					lines.add(line);
				}
			}
		}		

		checkSorting(lines);
	}

	/**
	 * Check that given lines are sorted correctly and throw and exception if that is not the case.
	 * 
	 * @param lines
	 * @throws UnsortedDataException
	 */
	private void checkSorting(List<String> lines) throws UnsortedDataException {
		
		Region previousRegion = null;
		
		for (String line : lines) {
			
			getParser().setLine(line);
			Region region = getParser().getRegion();
			
			if (previousRegion != null) {
				if (previousRegion.start.compareTo(region.start) > 0) {
					throw new UnsortedDataException("File " + getFile() + " isn't sorted correctly. " +
							"Please sort the file first.");
				}
			}
			previousRegion = region;
		}	
	}

	/**
	 * Initialize index with the location of first and last row.
	 * 
	 * @throws IOException
	 * @throws GBrowserException
	 */
	private void readEnds() throws IOException, GBrowserException {
		
		long firstLinePosition = 0;
		Region region = null;
		
		//Find the first line of real content
		do {
			//Add first row to index
			getFile().setLineReaderPosition(firstLinePosition);
			String firstLine = getFile().getNextLine();		
			getParser().setLine(firstLine);
			region = getParser().getRegion();
			
			if (region != null) {
				index.put(region.start, firstLinePosition);
			} else {
				firstLinePosition += firstLine.length() + 1; //plus one for new line character
			}
			
		} while (region == null);
		
		//Add last row to index
		String lastLine = getFile().getLastLine();		
		getParser().setLine(lastLine);
		region = getParser().getRegion();
		index.put(region.start, getFile().length() - lastLine.length() - 1);
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
									
			lines.add(line);
		}
		
		return lines;
	}
	
	public class LineIterator implements Iterator<String> {			
		
		private String next;

		public LineIterator() throws IOException, GBrowserException {
			BinarySearchIndex.this.getFile().setLineReaderPosition(0);
		}
		
		public void getNext() {
			if (next == null) {
				try {
					next = BinarySearchIndex.this.getFile().getNextLine();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		}

		@Override
		public boolean hasNext() {
			getNext();
			return next != null;
		}

		@Override
		public String next() {
			if (hasNext()) {
				String returnValue = next;
				next = null;
				return returnValue;				
			} else {
				throw new NoSuchElementException();
			}
		}

		@Override
		public void remove() {
			throw new NotImplementedException();
		}		
	}
	
	@Override
	public Iterator<String> getFileLineIterator() throws IOException, GBrowserException {
		
		return new LineIterator();	
	}

	public TreeMap<IndexKey, String> getFileLines(Region request) throws IOException, GBrowserException {
		
		if (request.start.compareTo(request.end) > 0) {
			//Illegal request, skip it
			//throw new IllegalArgumentException();
		}
		
		//Search start position
		
		//Decrease request start position by one to make sure that we get all lines when there  
		//are equal lines (in sort order) with the index entry in the file also before that index entry.
		BpCoord requestFloor = new BpCoord(request.start.bp - 1, request.start.chr);
		long floorFilePosition = binarySearch(requestFloor);
		
		getFile().setLineReaderPosition(floorFilePosition);
		
		TreeMap<IndexKey, String> lines = new TreeMap<IndexKey, String>();
		
		String line = null;
		long lineBytePosition = floorFilePosition;
		
		while ((line = getFile().getNextLine()) != null) {
			
			if ("".equals(line)) {
				//First byte was new line character
				continue;
			}
			
			getParser().setLine(line);
			Region region = getParser().getRegion();
			
			if (request.contains(region.start)) {
				lines.put(new IndexKey(region.start, lineBytePosition), line);
			}
			
			if (request.end.compareTo(region.start) < 0) {
				break;
			}
			
			lineBytePosition += line.length() + 1; // plus one for new line character
		}
		
		return lines;
	}
	
	/**
	 * Search previous and following index entries around the requested position from the index. Create a new 
	 * index entry between these entries, if the distance between them is too big for line-by-line reading.
	 * Effectively makes a binary search on the file to find a location that is reasonably close to request position
	 * and precedes it.  
	 * 
	 * @param position
	 * @return
	 * @throws IOException
	 * @throws GBrowserException
	 */
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

	/**
	 * Create a new index entry between the given file positions.
	 * 
	 * @param floorFilePosition
	 * @param ceilingFilePosition
	 * @throws IOException
	 * @throws GBrowserException
	 */
	private void splitIndexRegion(long floorFilePosition,
			long ceilingFilePosition) throws IOException, GBrowserException {
		
		long centerFilePosition = (floorFilePosition + ceilingFilePosition) / 2;
		
		getFile().setLineReaderPosition(centerFilePosition);
		//Skip the first line, because we don't know if it is complete or partial 		
		long partialLineLength = getFile().getNextLine().length() + 1; //plus one because of new line character
		String line = getFile().getNextLine();		
		
		getParser().setLine(line);
		Region region = getParser().getRegion();

		//index file positions point to preceding new line character
		index.put(region.start, centerFilePosition + partialLineLength + line.length());
	}

	public RandomAccessLineDataSource getFile() {
		return file;
	}

	public LineParser getParser() {
		return parser;
	}
}
