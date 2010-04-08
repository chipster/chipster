package fi.csc.microarray.databeans;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
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
import fi.csc.microarray.databeans.handlers.DataBeanHandler;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.InputStreamSource;
import fi.csc.microarray.util.StreamStartCache;

/**
 * <p>DataBean is the basic unit of databeans package. It holds a chunk
 * of content (like file) and related metadata. A good way to think
 * about databeans is to view them as a high level substitute for files.
 * DataBeans are managed by a DataManager: every bean has exactly one
 * manager. Beans should not be mixed across manager boundaries.</p>
 * 
 * <p>DataBean stores two kinds of relationships: hierarchical relationships 
 * of DataBeans and DataFolders and link relationships of related beans.</p>
 * 
 * @author Aleksi Kallio
 *
 */
public class DataBean extends DataItemBase {
	
	/**
	 * Traversal specifies the way of traversing links. 
	 *
	 */
	public enum Traversal {
		DIRECT,
		REVERSED,
		BIDIRECTIONAL;
		
		public boolean isDirect() {
			return this == DIRECT || this == BIDIRECTIONAL;
		}
		
		public boolean isReversed() {
			return this == REVERSED || this == BIDIRECTIONAL;
		}
	}
	
	/**
	 * Link represents a relationship between two beans.
	 *
	 */
	public enum Link {
		/**
		 * Relationship where one bean describes (acts as metadata for) another.
		 */
		ANNOTATION("Annotation"),
		
		/**
		 * Relationship where other bean has been used to derive another (i.e. beans are an input and an output of an operation).
		 */
		DERIVATION("Derivation"),
		
		/**
		 * Relationship where other bean has been modified into another (i.e. user has modified and saved a new copy of a bean).
		 * Philosophically modification is a special case of derivation, but they are
		 * modelled as two different types to separate manual (modification) and 
		 * automatic production (derivation) of beans.  
		 */
		MODIFICATION("Modification"),
		
		/**
		 * Relationship where two beans belong to a same group.
		 */
		GROUPING("Grouping");
		
		private String name;
		
		private Link(String name) {
			this.name = name;
		}
		
		public String toString(){
			return this.name;
		}

		public static Link[] userEditableValues() {
			return new Link[] {ANNOTATION, DERIVATION, MODIFICATION};
		}

		public static Link[] derivationalTypes() {
			return new Link[] {DERIVATION, MODIFICATION};
		}
	}
	
	
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
	private static final Logger logger = Logger.getLogger(DataBean.class);

	protected DataManager dataManager;
	protected StreamStartCache streamStartCache = null;
	private HashMap<String, Object> contentCache = new HashMap<String, Object>();
	
	private URL cacheUrl = null;
	private boolean contentChanged = true;
	private Lock contentLock = new ReentrantLock();

	private LinkedList<LinkedBean> outgoingLinks = new LinkedList<LinkedBean>();
	private LinkedList<LinkedBean> incomingLinks = new LinkedList<LinkedBean>();
	
	protected Date date;
	private Operation sourceOperation;
	private String notes;

	protected ContentType contentType;

	
	private URL contentUrl;
	private DataBeanHandler handler;
	

	private File contentFile;

	
	
	public DataBean(String name, ContentType contentType, Date date, DataBean[] sources, DataFolder parentFolder, DataManager manager, File contentFile) {
		
		this.dataManager = manager;
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


	public DataBean(String name, URL contentUrl, ContentType contentType, Date date, DataBean[] sources, DataFolder parentFolder, DataManager manager, DataBeanHandler handler) {
		
		this.dataManager = manager;
		this.name = name;
		this.contentUrl = contentUrl;
		this.handler = handler;
		this.date = date;
		this.parent = parentFolder;
		
		// add this as parent folders child
		if (parentFolder != null) {
			parentFolder.addChild(this);
		}
		
		for (DataBean source : sources) {
			source.addLink(Link.DERIVATION, this);
		}

		this.contentType = contentType;
	}

	
	
	
	
	
	

	
	
	
	
	
	
	
	
	


	




	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	





	/**
	 * @return A String presentation of this dataset (in this case,
	 * 		   simply the name, to be shown on e.g. the tree).
	 */
	public String toString() {
		return getName();
	}
	
	
	/**
	 * Checks if bean's content type is any of the given alternatives. Case insensitive.
	 */
	public boolean isContentTypeCompatitible(String... alternatives) {
		for (String alternative : alternatives) {
			if (alternative.toLowerCase().equals(getContentType().getType().toLowerCase())) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Returns result of the feature query. Acts as a gateway to Feature API which
	 * is the way bean contents can be used besides using the raw binary content.
	 * Available features depend on what feature factories are plugged to DataManager.
	 */
	public QueryResult queryFeatures(String request) {
		Feature feature = new RequestExecuter(dataManager).execute(request, this);
		if (feature == null) {
			throw new UnsupportedOperationException("request " + request + " not possible from " + this.getName());
		}
		return new QueryResult(feature);
	}
	
	
	/**
	 * Returns raw binary content of this bean. Use Feature API for 
	 * higher level accessing.
	 * 
	 * @see #queryFeatures(String)
	 */
	public InputStream getContentByteStream() throws IOException {
		if (streamStartCache != null) {
			return streamStartCache.getInputStream();
		} else {
			logger.debug("using non-cached stream");
			return getRawContentByteStream();
		}
	}
	
	/**
	 *  TODO should be integrated and hidden away
	 * @throws MicroarrayException
	 * @throws IOException
	 */
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
	
	
	
	/**
	 * Convenience method which includes all beans into result.
	 * 
	 * @see #traverseLinks(fi.csc.microarray.databeans.DataBean.Link, fi.csc.microarray.databeans.DataBean.Traversal, DataBeanSelector)
	 */
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

	
	
	/**
	 * Does a breadth first search following links of given type. Uses selector
	 * to decide if beans should be included into result and if they should be traversed
	 * further. Selected beans are returned in the order (breadth first) they were encountered.
	 * No duplicate beans are returned.
	 * 
	 * @param type
	 * @param traversal
	 * @param selector
	 * @return
	 */
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
	
	
	/**
	 * Puts named object to content cache. Content cache is a hashmap type per DataBean cache where 
	 * objects can be stored. Cache is emptied every time bean content is changed, so it is suited
	 * for caching results that are derived from contents of a single bean. Cache is not persistant,
	 * and generally, user should never assume cached values to be found.  
	 */
	public void putToContentCache(String name, Object value) {
		this.contentCache.put(name, value);
	}
	
	/**
	 * Gets named object from content cache.
	 * 
	 * @see #addToContentCache(String, Object)
	 * @param name
	 * @return
	 */
	public Object getFromContentCache(String name) {
		return this.contentCache.get(name);
	}
	
	protected void resetContentCache() {
		this.contentCache.clear();
	}
	
	/**
	 * Acquire the lock which guards the content from being changed.
	 * 
	 * Needed for example when transferring contents to file broker.
	 * 
	 */
	public void lockContent() {
		this.contentLock.lock();
	}
	
	/**
	 * Release the content lock.
	 * 
	 */
	public void unlockContent() {
		this.contentLock.unlock();
	}

	/**
	 * Indicate whether the contents have been changed since the contents
	 * were last time uploaded to file broker.
	 * 
	 */
	public boolean isContentChanged() {
		return this.contentChanged;
	}

	
	/**
	 * Set content changed status. Shoulc be called with true everytime
	 * content is changed.
	 * 
	 * @param contentChanged
	 */
	public void setContentChanged(boolean contentChanged) {
		this.contentChanged = contentChanged;
	}
	
	/**
	 * Get the location of the remote copy for the content file. 
	 * Usually the copy is located at the file broker.
	 * 
	 * @return may be null
	 */
	public URL getCacheUrl() {
		return this.cacheUrl;
	}
	
	/**
	 * Set the location of the remote copy of the content file.
	 * Usually the copy is located at the file broker.
	 * 
	 * @param url
	 */
	public void setCacheUrl(URL url) {
		this.cacheUrl = url;
	}

	/**
	 * @return The date and time when this dataset was created.
	 */
	public Date getDate() {
		return this.date;
	}
	
	public void setName(String newName) {
		super.setName(newName);

		ContentChangedEvent cce = new ContentChangedEvent(this);		
		dataManager.dispatchEventIfVisible(cce);
	}

	
	/**
	 * Removes a link. Does not remove links of different type between same beans.
	 */
	public void removeLink(Link type, DataBean target) {
		for (LinkedBean outgoingLink : outgoingLinks) {
			if (outgoingLink.bean == target && outgoingLink.link == type) {
				
				outgoingLinks.remove(outgoingLink);
				
				DataBean bean = outgoingLink.bean;
				
				for (LinkedBean incomingLink : bean.incomingLinks) {
					if (incomingLink.bean == this && incomingLink.link == type) {
						bean.incomingLinks.remove(incomingLink);
						
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
	

	/**
	 * Creates a link between this and target bean. Links represent relationships
	 * between beans. Links have hardcoded types.
	 * 
	 * @see Link
	 */
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
		target.incomingLinks.add(new LinkedBean(type, this));
		outgoingLinks.add(new LinkedBean(type, target));

		// dispatch events only if both visible
		if (this.getParent() != null && target.getParent() != null) {
			LinksChangedEvent lce = new LinksChangedEvent(this, target, type, true);
			dataManager.dispatchEvent(lce);
		}

	}
	
	
	/**
	 * @return outgoing links of given type.
	 */
	public List<DataBean> getLinkTargets(Link... types) {
		return getLinkedBeans(types, outgoingLinks);
	}

	/**
	 * @return incoming links of given type.
	 */
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

	/**
	 * @return MIME content type of the bean.
	 */
	public ContentType getContentType() {
		return contentType;
	}

	public void setContentType(ContentType contentType) {
		this.contentType = contentType;
	}

	
	public InputStream getRawContentByteStream() throws IOException {
		InputStream is;
		try {
			if (this.contentUrl != null) {
				is = handler.getInputStream(this);
			} else {
				is = new FileInputStream(contentFile);
			}
		} catch (IOException e) {
			throw e;
		}
		return is;
	}

	/**
	 * A convenience method for gathering streamed binary content into
	 * a byte array.
	 * 
	 *   @see #getContentByteStream()
	 */
	public byte[] getContents() throws MicroarrayException {
		try {
			return Files.inputStreamToBytes(this.getContentByteStream());
		} catch (IOException e) {
			throw new MicroarrayException(e);
		}
		
	}

	/**
	 * 
	 * 	
	 * Allows rewriting of raw bean content. Close the stream by calling closeContentOutputStream(...)
	 * on the same bean.
	 * @see #closeContentOutputStreamAndUnlockDataBean(OutputStream)
	 *
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


	/**
	 * Closes output stream and generates required events.
	 */
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
	
	/**
	 * Returns content size in bytes.
	 */
	public long getContentLength() {
		if (this.contentUrl != null) {
			try {
				return handler.getContentLength(this);
			} catch (IOException e) {
				// FIXME
				throw new RuntimeException(e);
			}
		} else {
			return contentFile.length();
		}
	}


	public URL getContentUrl() {
		return contentUrl;
	}


	public void setContentUrl(URL contentUrl) {
		this.contentUrl = contentUrl;
	}

	/**
	 * Should only be used when saving the bean inside a session.
	 * 
	 * @param handler
	 */
	public void setHandler(DataBeanHandler handler) {
		this.handler = handler;
	}







}