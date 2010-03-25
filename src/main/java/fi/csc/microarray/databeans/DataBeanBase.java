package fi.csc.microarray.databeans;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.QueryResult;
import fi.csc.microarray.databeans.features.RequestExecuter;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.InputStreamSource;
import fi.csc.microarray.util.StreamStartCache;
import fi.csc.microarray.util.Strings;

public abstract class DataBeanBase extends DataItemBase implements DataBean {
	
	private static class LinkedBean {
		
		LinkedBean(Link link, DataBean bean) {
			this.link = link;
			this.bean = bean;			
		}
		
		Link link;
		DataBean bean;
	}
	
	
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(DataBeanBase.class);

	protected DataManagerBase dataManager;
	protected StreamStartCache streamStartCache = null;
	private HashMap<String, Object> contentCache = new HashMap<String, Object>();
	
	private URL url = null;
	private boolean contentChanged = true;
	private Lock contentLock = new ReentrantLock();

	private LinkedList<LinkedBean> outgoingLinks = new LinkedList<LinkedBean>();
	private LinkedList<LinkedBean> incomingLinks = new LinkedList<LinkedBean>();
	
	protected Date date;
	private Operation sourceOperation;
	private String notes;

	protected ContentType contentType;

	
	protected DataBeanBase(DataManagerBase dataManager) {
		this.dataManager = dataManager;
	}

	public abstract InputStream getRawContentByteStream() throws MicroarrayException;	

	/**
	 * @return A String presentation of this dataset (in this case,
	 * 		   simply the name, to be shown on e.g. the tree).
	 */
	public String toString() {
		return getName();
	}
	
	public boolean isContentTypeCompatitible(String... alternatives) {
		for (String alternative : alternatives) {
			if (alternative.toLowerCase().equals(getContentType().getType().toLowerCase())) {
				return true;
			}
		}
		return false;
	}
	
	public QueryResult queryFeatures(String request) {
		Feature feature = new RequestExecuter(dataManager).execute(request, this);
		if (feature == null) {
			throw new UnsupportedOperationException("request " + request + " not possible from " + this.getName());
		}
		return new QueryResult(feature);
	}
	
	public InputStream getContentByteStream() throws MicroarrayException, IOException {
		if (streamStartCache != null) {
			return streamStartCache.getInputStream();
		} else {
			logger.debug("using non-cached stream");
			return getRawContentByteStream();
		}
	}
	
	public void initialiseStreamStartCache() throws MicroarrayException, IOException {

		this.streamStartCache = new StreamStartCache(getRawContentByteStream(), new InputStreamSource() {
			public InputStream getInputStream() {
				try {
					return getRawContentByteStream();
				} catch (Exception e) {
					throw new RuntimeException(e);
				}
			}
		});
	}
	
	public List<DataBean> traverseLinks(Link[] types, Traversal traversal) {
		
		DataBeanSelector acceptAllSelector = new DataBeanSelector() {
			public boolean shouldSelect(DataBean bean) {
				return true;
			}
			public boolean shouldTraverse(DataBean bean) {
				return true;
			}
		};
		
		return traverseLinks(types, traversal, acceptAllSelector);
	}

	public List<DataBean> traverseLinks(Link[] types, Traversal traversal, DataBeanSelector selector) {
		LinkedList<DataBean> selected = new LinkedList<DataBean>();
		LinkedList<DataBean> traversed = new LinkedList<DataBean>();
		traversed.add(this);
		conditionallySelect(selector, selected, this);
		
		LinkedHashSet<DataBean> newBeans;
		do {
			newBeans = new LinkedHashSet<DataBean>();
		
			for (DataBean bean : traversed) {
				LinkedHashSet<DataBean> linkedBeans = new LinkedHashSet<DataBean>();
				
				if (traversal.isDirect()) {
					for (Link type : types) {
						linkedBeans.addAll(bean.getLinkTargets(type));
					}
				}
				
				if (traversal.isReversed()) {
					for (Link type : types) {
						linkedBeans.addAll(bean.getLinkSources(type));
					}
				}
			
				for (DataBean linkedBean : linkedBeans) {
					if (!traversed.contains(linkedBean) && selector.shouldTraverse(linkedBean)) {
						newBeans.add(linkedBean);
					}
				}
			}
			
			for (DataBean newBean : newBeans) {
				traversed.add(newBean);
				conditionallySelect(selector, selected, newBean);

			}
			
		} while (newBeans.size() > 0); // iterate as long as we make progress
		
		return selected;
	}

	private void conditionallySelect(DataBeanSelector selector, LinkedList<DataBean> selected, DataBean bean) {
		if (!selected.contains(bean) && selector.shouldSelect(bean)) {
			selected.add(bean);
		}
	}
	
	public void putToContentCache(String name, Object value) {
		this.contentCache.put(name, value);
	}
	
	public Object getFromContentCache(String name) {
		return this.contentCache.get(name);
	}
	
	protected void resetContentCache() {
		this.contentCache.clear();
	}
	
	public void lockContent() {
		this.contentLock.lock();
	}
	
	public void unlockContent() {
		this.contentLock.unlock();
	}

	public boolean isContentChanged() {
		return this.contentChanged;
	}

	public void setContentChanged(boolean contentChanged) {
		this.contentChanged = contentChanged;
	}
	
	public URL getUrl() {
		return this.url;
	}
	public void setUrl(URL url) {
		this.url = url;
	}

	public Date getDate() {
		return this.date;
	}
	
	public void setName(String newName) {
		super.setName(newName);

		ContentChangedEvent cce = new ContentChangedEvent(this);		
		dataManager.dispatchEventIfVisible(cce);
	}

	public void removeLink(Link type, DataBean target) {
		for (LinkedBean outgoingLink : outgoingLinks) {
			if (outgoingLink.bean == target && outgoingLink.link == type) {
				
				outgoingLinks.remove(outgoingLink);
				
				// FIXME don't reference incoming directly
				DataBeanBase fsBean = ((DataBeanBase)outgoingLink.bean);
				
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
		// FIXME don't reference incomingLinks directly
		((DataBeanBase)target).incomingLinks.add(new LinkedBean(type, this));
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

	
	public void setNotes(String notes) {
		this.notes = notes;		
	}

	public String getNotes() {
		return notes;
	}

	public ContentType getContentType() {
		return contentType;
	}

	public void setContentType(ContentType contentType) {
		this.contentType = contentType;
	}

	
	
}
