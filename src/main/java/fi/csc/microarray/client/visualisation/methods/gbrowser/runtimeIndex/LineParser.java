package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;


/**
 * Parsers for different file formats have to implement these methods to offer an unified
 * interface for Indexes. 
 * 
 * @author klemela
 */
public interface LineParser {

	public Region getRegion();

	public boolean setLine(String line);

	public boolean isContentLine();

	public FileLine getFileLine();	
}
