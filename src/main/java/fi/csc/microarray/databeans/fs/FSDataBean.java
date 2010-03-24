package fi.csc.microarray.databeans.fs;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Date;

import fi.csc.microarray.databeans.ContentChangedEvent;
import fi.csc.microarray.databeans.ContentType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBeanBase;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataManagerBase;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.Files;

public class FSDataBean extends DataBeanBase {
	
	
	private File contentFile;
	
	
	public FSDataBean(String name, ContentType contentType, Date date, DataBean[] sources, DataFolder parentFolder, DataManagerBase manager, File contentFile) {
		super(manager);

		this.name = name;
		this.date = date;
		this.parent = parentFolder;
		
		// add this as parent folders child
		if (parentFolder != null) {
			parentFolder.addChild(this);
		}
		
		for (DataBean source : sources) {
			source.addLink(Link.DERIVATION, this);
		}

		this.contentFile = contentFile;
		this.contentType = contentType;
	}
		
	
	
	public InputStream getRawContentByteStream() throws MicroarrayException {
		InputStream is;
		try {
			is = new FileInputStream(this.contentFile);
		} catch (FileNotFoundException e) {
			throw new MicroarrayException(e);
		}
		return is;
	}

	public byte[] getContents() throws MicroarrayException {
		try {
			return Files.inputStreamToBytes(this.getContentByteStream());
		} catch (IOException e) {
			throw new MicroarrayException(e);
		}
		
	}

	/**
	 * Returns OutputStream that can be used to rewrite this bean's contents. 
	 * Calling this method results in disabling caching for this bean.
	 */
	public OutputStream getContentOutputStreamAndLockDataBean() throws MicroarrayException, IOException {
		lockContent();
		setContentChanged(true);
		this.streamStartCache = null; // caching is disabled
		resetContentCache();
		OutputStream os;
		try {
			os = new FileOutputStream(this.contentFile);
		} catch (FileNotFoundException e) {
			throw new MicroarrayException(e);
		}
		return os;
	}

	public void closeContentOutputStreamAndUnlockDataBean(OutputStream out)
			throws MicroarrayException, IOException {
		
		try {
			out.close();
		} finally {
			unlockContent();
		}
		ContentChangedEvent cce = new ContentChangedEvent(this);
		dataManager.dispatchEventIfVisible(cce);
	}

	public void delete() {
		lockContent();
		try {
			this.contentFile.delete();
			this.contentFile = null;
			this.contentType = null;			
		} finally {
			unlockContent();
		}
	}
	
	public long getContentLength() {
		return contentFile.length();
	}

	public File getContentFile() {
		return contentFile;
	}
}
