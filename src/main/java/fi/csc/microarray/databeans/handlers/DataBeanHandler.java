package fi.csc.microarray.databeans.handlers;

import java.io.IOException;
import java.io.InputStream;

import fi.csc.microarray.databeans.DataBean;

public interface DataBeanHandler {

	public InputStream getInputStream(DataBean dataBean) throws IOException;
	
	public long getContentLength(DataBean dataBean) throws IOException;
	
	public void delete(DataBean dataBean);
}
