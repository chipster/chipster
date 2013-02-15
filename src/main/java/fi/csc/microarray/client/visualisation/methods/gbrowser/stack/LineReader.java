package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.IOException;

public interface LineReader {

	public void setPosition(long position) throws IOException;

	public String readLine() throws IOException;

	public void close();
	
	public long length() throws IOException;
}
