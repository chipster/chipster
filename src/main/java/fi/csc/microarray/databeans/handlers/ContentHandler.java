package fi.csc.microarray.databeans.handlers;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import fi.csc.microarray.databeans.DataBean;

public interface ContentHandler {

	public InputStream getInputStream(DataBean dataBean) throws IOException;
	
	public OutputStream getOutputStream(DataBean dataBean) throws IOException;
	
	public long getContentLength(DataBean dataBean) throws IOException;
	
	public void delete(DataBean dataBean);
}
