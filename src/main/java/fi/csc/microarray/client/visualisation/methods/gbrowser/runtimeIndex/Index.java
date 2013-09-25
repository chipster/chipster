package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.IOException;
import java.util.Iterator;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.IndexKey;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

public abstract class Index {

	/**
	 * Get lines from file. Only lines with start position within the requestRegion are returned, request start position
	 * is inclusive, end position is exclusive.
	 * 
	 * @param requestRegion
	 * @return
	 * @throws IOException
	 * @throws GBrowserException
	 */
	public abstract TreeMap<IndexKey, String> getFileLines(Region requestRegion) throws IOException, GBrowserException;

	public abstract Iterator<String> getFileLineIterator() throws IOException, GBrowserException;
}
