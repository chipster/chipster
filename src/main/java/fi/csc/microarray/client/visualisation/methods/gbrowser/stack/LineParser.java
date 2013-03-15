package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;


public interface LineParser {

	public Region getRegion();

	public boolean setLine(String line);

	public boolean isContentLine();
}
