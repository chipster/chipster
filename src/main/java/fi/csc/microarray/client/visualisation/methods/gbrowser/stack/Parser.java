package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;


public interface Parser {

	public Region getRegion();

	public boolean setLine(String line);
}
