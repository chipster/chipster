package fi.csc.microarray.client.visualisation.methods.gbrowser.util;

import org.broad.tribble.readers.TabixReader;
import org.broad.tribble.readers.TabixReader.Iterator;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class TabixUtil {
	public static Iterator getTabixIterator(TabixDataSource dataSource, Region request) {
		String chromosome = request.start.chr.toNormalisedString();

		//limit to integer range
		int start = (int) Math.min(Integer.MAX_VALUE, request.start.bp);
		int end = (int) Math.min(Integer.MAX_VALUE, request.end.bp);
		
		//Extend area to be able to draw introns at screen edge, but don't go over MAX_VALUE, or below 1
		//TODO Be more clever to avoid getting so much useless data
		int EXTRA = 500000; //O,5M should be enought for the longest human introns http://www.bioinfo.de/isb/2004040032/
		
		start = (int) Math.max((long)start - EXTRA, 1);
		end = (int) Math.min((long)end + EXTRA, Integer.MAX_VALUE);

		//Check that region is below max bin size of Tabix
		int MAX_BIN_SIZE = 512*1024*1024 - 2;

		start = (int) Math.min(MAX_BIN_SIZE, start);
		end = (int) Math.min(MAX_BIN_SIZE, end);

		start = (int) Math.max(1, start);
		end = (int) Math.max(1, end);

		String queryRegion = chromosome + ":" + start + "-" + end;

		TabixReader.Iterator iter = null;
		
		try {
			iter = dataSource.getReader().query(queryRegion);
		} catch (ArrayIndexOutOfBoundsException e) {
			//No such chromosome
		}
			
		return iter;
	}
}
