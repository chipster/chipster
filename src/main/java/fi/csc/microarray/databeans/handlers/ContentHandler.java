package fi.csc.microarray.databeans.handlers;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import fi.csc.microarray.databeans.DataBean.ContentLocation;

public interface ContentHandler {

	public InputStream getInputStream(ContentLocation location) throws IOException;
	
	public OutputStream getOutputStream(ContentLocation location) throws IOException;
	
	public long getContentLength(ContentLocation location) throws IOException;
	
	public void canBeDeleted(ContentLocation location);
}
