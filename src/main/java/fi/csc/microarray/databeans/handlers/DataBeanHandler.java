package fi.csc.microarray.databeans.handlers;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import fi.csc.microarray.databeans.Dataset;

public interface DataBeanHandler {

	public InputStream getInputStream(Dataset dataBean) throws IOException;
	
	public OutputStream getOutputStream(Dataset dataBean) throws IOException;
	
	public long getContentLength(Dataset dataBean) throws IOException;
	
	public void delete(Dataset dataBean);
}
