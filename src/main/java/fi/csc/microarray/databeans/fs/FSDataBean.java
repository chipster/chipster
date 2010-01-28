package fi.csc.microarray.databeans.fs;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Date;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.databeans.ContentChangedEvent;
import fi.csc.microarray.databeans.ContentType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBeanBase;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataManagerBase;
import fi.csc.microarray.databeans.LinksChangedEvent;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.Files;

public class FSDataBean extends DataBeanBase {
	
	private static class LinkedBean {
		
		LinkedBean(Link link, DataBean bean) {
			this.link = link;
			this.bean = bean;			
		}
		
		Link link;
		DataBean bean;
	}
	
	private ContentType contentType;
	private String name;
	private Operation sourceOperation;
	private Date date;
	private String notes;
	
	private DataFolder parentFolder;
	private LinkedList<LinkedBean> outgoingLinks = new LinkedList<LinkedBean>();
	private LinkedList<LinkedBean> incomingLinks = new LinkedList<LinkedBean>();

	private File contentFile;
	
	
	public FSDataBean(String name, ContentType contentType, Date date, DataBean[] sources, DataFolder parentFolder, DataManagerBase manager, File contentFile) {
		super(manager);

		this.name = name;
		this.date = date;
		this.parentFolder = parentFolder;
		
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

	public Date getDate() {
		return this.date;
	}

	public String getName() {
		return this.name;
	}
	
	public void removeLink(Link type, DataBean target) {
		for (LinkedBean outgoingLink : outgoingLinks) {
			if (outgoingLink.bean == target && outgoingLink.link == type) {
				
				outgoingLinks.remove(outgoingLink);
				FSDataBean fsBean = ((FSDataBean)outgoingLink.bean);
				
				for (LinkedBean incomingLink : fsBean.incomingLinks) {
					if (incomingLink.bean == this && incomingLink.link == type) {
						fsBean.incomingLinks.remove(incomingLink);
						
						// both links were found
						LinksChangedEvent lce = new LinksChangedEvent(this, target, type, false);
						dataManager.dispatchEventIfVisible(lce);
						return; 
					}
				}
			}
		}
		throw new RuntimeException("internal error: failed locate links for: " + this.getName()  + " <" + type + "> " + target.getName());
	}
	
	public DataFolder getParent() {
		return this.parentFolder;
	}

	public void setName(String newName) {
		this.name = newName;

		ContentChangedEvent cce = new ContentChangedEvent(this);		
		dataManager.dispatchEventIfVisible(cce);
	}

	public void addLink(Link type, DataBean target) {
		if (target == null) {
			return;
		}		
		
		// FIXME add more internal state validation to FSDataBean
//		for (LinkedBean linkedBean : outgoingLinks) {
//			if (linkedBean.bean == target && linkedBean.link == type) {
//				throw new RuntimeException("duplicate link");
//			}
//		}

		// make both parties aware of the link
		((FSDataBean)target).incomingLinks.add(new LinkedBean(type, this));
		outgoingLinks.add(new LinkedBean(type, target));

		// dispatch events only if both visible
		if (this.getParent() != null && target.getParent() != null) {
			LinksChangedEvent lce = new LinksChangedEvent(this, target, type, true);
			dataManager.dispatchEvent(lce);
		}

	}
	
		
	/**
	 * @return The parent dataset of this dataset.
	 */
	public List<DataBean> getLinkTargets(Link... types) {
		return getLinkedBeans(types, outgoingLinks);
	}

	public List<DataBean> getLinkSources(Link... types) {
		return getLinkedBeans(types, incomingLinks);
	}
	
	private List<DataBean> getLinkedBeans(Link[] types, LinkedList<LinkedBean> links) {
		LinkedList<DataBean> targets = new LinkedList<DataBean>();
		for (LinkedBean linkTarget : links) {
			for (Link type : types) {
				if (linkTarget.link == type) {
					targets.add(linkTarget.bean);
					break;
				}
			}
		}
		return targets;
	}
	
	/**
	 * @return The operation that has been selected for this dataset (and which
	 * 		   may have already been done at least once to produce another
	 * 		   dataset, or which has not yet been conducted, and maybe never
	 * 		   will, depending on the user).
	 */
	public Operation getOperation() {
		return sourceOperation;
	}

	/**
	 * Associates the given operation with this DataBean.
	 * 
	 * @param operation to associate
	 */
	public void setOperation(Operation operation) {
		this.sourceOperation = operation;
	}
	
	public ContentType getContentType() {
		return contentType;
	}
	
	public void setContentType(ContentType contentType) {
		this.contentType = contentType;
	}

	public void setNotes(String notes) {
		this.notes = notes;		
	}

	public String getNotes() {
		return notes;
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

	public void setParent(FSDataFolder dataFolder) {
		this.parentFolder = dataFolder;		
	}
	
}
