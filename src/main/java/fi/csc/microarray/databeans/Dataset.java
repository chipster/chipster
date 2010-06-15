package fi.csc.microarray.databeans;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.locks.ReentrantReadWriteLock;

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
public class Dataset extends DataItemBase {
	
	public enum DataBeanType {
		LOCAL_USER,
		LOCAL_TEMP,
		LOCAL_SESSION;
	}
	
	
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
		
		LinkedBean(Link link, Dataset bean) {
			this.link = link;
			this.bean = bean;			
		}
		
		Link link;
		Dataset bean;
	}
	
	
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(Dataset.class);

	protected DataManager dataManager;
	protected StreamStartCache streamStartCache = null;
	private HashMap<String, Object> contentCache = new HashMap<String, Object>();
	
	private URL cacheUrl = null;
	private boolean contentChanged = true;
	private ReentrantReadWriteLock lock = new ReentrantReadWriteLock(true);


	private LinkedList<LinkedBean> outgoingLinks = new LinkedList<LinkedBean>();
	private LinkedList<LinkedBean> incomingLinks = new LinkedList<LinkedBean>();
	
	protected Date date;
	private Operation sourceOperation;
	private String notes;

	protected ContentType contentType;

	
	private DataBeanType type;
	private String repositoryName;
	private URL url;
	private DataBeanHandler handler;


	public Dataset(String name, DataBeanType type, String repositoryName, URL contentUrl, ContentType contentType, Date date, Dataset[] sources, DataFolder parentFolder, DataManager manager, DataBeanHandler handler) {
		
		this.dataManager = manager;
		this.name = name;
		this.url = contentUrl;
		this.type = type;
		this.repositoryName = repositoryName;
		this.handler = handler;
		this.date = date;
		this.parent = parentFolder;
		
		
		// add this as parent folders child
		if (parentFolder != null) {
			parentFolder.addChild(this);
		}
		
		for (Dataset source : sources) {
			source.addLink(Link.DERIVATION, this);
		}

		this.contentType = contentType;
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
	 * Returns result of the feature query. Acts as a gateway to Feature API which
	 * is the way bean contents can be used besides using the raw binary content.
	 * Available features depend on what feature factories are plugged to DataManager.
	 */
	public QueryResult queryFeatures(String request) {
		try {
			lock.readLock().lock();
			Feature feature = new RequestExecuter(dataManager).execute(request, this);
			if (feature == null) {
				throw new UnsupportedOperationException("request " + request + " not possible from " + this.getName());
			}
			return new QueryResult(feature);
		} finally {
			lock.readLock().unlock();
		}
	}



	/**
	 * Returns raw binary content of this bean. Use Feature API for 
	 * higher level accessing.
	 * 
	 * @see #queryFeatures(String)
	 * 
	 * FIXME change name, locks stream start cache
	 */
	public InputStream getContentByteStream() throws IOException {
		return getRawContentByteStream();

		//lock.readLock().lock();
//		if (streamStartCache != null) {
//			return streamStartCache.getInputStream();
//		} else {
//			logger.debug("using non-cached stream");
//			return getRawContentByteStream();
//		}
	}



	/**
	 * A convenience method for gathering streamed binary content into
	 * a byte array.
	 * @throws IOException 
	 * 
	 *   @see #getContentByteStream()
	 */
	public byte[] getContents() throws IOException {
		try {
			return Files.inputStreamToBytes(this.getContentByteStream());
		} finally {
			// FIXME release inputstream
		}
		
	}
	
	public byte[] getContents(long maxLength) throws IOException {
		ByteArrayOutputStream outputStream = new ByteArrayOutputStream();
		InputStream in = this.getContentByteStream();
		long readCount = 0;
		
		for (int b = in.read(); b != -1 && readCount < maxLength; b = in.read()) {
			outputStream.write(b);
			readCount++;
		}
		
		return outputStream.toByteArray();
	}



	/**
	 * Returns content size in bytes.
	 */
	public long getContentLength() {
		if (this.handler == null) {
			throw new IllegalStateException("Handler is null.");
		}
	
		try {
			return handler.getContentLength(this);
		} catch (IOException e) {
			// FIXME what to do?
			throw new RuntimeException(e);
		}
	}



	public void delete() {
//		lock.writeLock().lock();
		try {
			this.handler.delete(this);
			this.contentType = null;			
		} finally {
//			lock.writeLock().unlock();
		}
	}



	/**
	 * Puts named object to content cache. Content cache is a hashmap type per DataBean cache where 
	 * objects can be stored. Cache is emptied every time bean content is changed, so it is suited
	 * for caching results that are derived from contents of a single bean. Cache is not persistent,
	 * and generally, user should never assume cached values to be found.  
	 *
	 *	FIXME think about locking
	 *
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
		try {
//			this.lock.readLock().lock();
	
			return this.contentCache.get(name);
		} finally {
//			this.lock.readLock().unlock();
		}
	}



	/** 
	 * FIXME Think about locking.
	 */
	protected void resetContentCache() {
		this.contentCache.clear();
	}



	/**
	 *  TODO should be integrated and hidden away
	 * @throws MicroarrayException
	 * @throws IOException
	 */
	public void initialiseStreamStartCache() throws MicroarrayException, IOException {
		try {
//			this.lock.readLock().lock();
			this.streamStartCache = new StreamStartCache(getRawContentByteStream(), new InputStreamSource() {
				public InputStream getInputStream() {
					try {
						return getRawContentByteStream();
					} catch (Exception e) {
						throw new RuntimeException(e);
					}
				}
			});
		} finally {
//			this.lock.readLock().unlock();
		}
	}
	
	
	
	/**
		 * Creates a link between this and target bean. Links represent relationships
		 * between beans. Links have hardcoded types.
		 * 
		 * @see Link
		 */
		public void addLink(Link type, Dataset target) {
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
	 * Removes a link. Does not remove links of different type between same beans.
	 */
	public void removeLink(Link type, Dataset target) {
		for (LinkedBean outgoingLink : outgoingLinks) {
			if (outgoingLink.bean == target && outgoingLink.link == type) {
				
				outgoingLinks.remove(outgoingLink);
				
				Dataset bean = outgoingLink.bean;
				
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
	 * Convenience method which includes all beans into result.
	 * 
	 * @see #traverseLinks(fi.csc.microarray.databeans.Dataset.Link, fi.csc.microarray.databeans.Dataset.Traversal, DataBeanSelector)
	 */
	public List<Dataset> traverseLinks(Link[] types, Traversal traversal) {
		
		DataBeanSelector acceptAllSelector = new DataBeanSelector() {
			public boolean shouldSelect(Dataset bean) {
				return true;
			}
			public boolean shouldTraverse(Dataset bean) {
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
	public List<Dataset> traverseLinks(Link[] types, Traversal traversal, DataBeanSelector selector) {
		LinkedList<Dataset> selected = new LinkedList<Dataset>();
		LinkedList<Dataset> traversed = new LinkedList<Dataset>();
		traversed.add(this);
		conditionallySelect(selector, selected, this);
		
		LinkedHashSet<Dataset> newBeans;
		do {
			newBeans = new LinkedHashSet<Dataset>();
		
			for (Dataset bean : traversed) {
				LinkedHashSet<Dataset> linkedBeans = new LinkedHashSet<Dataset>();
				
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
			
				for (Dataset linkedBean : linkedBeans) {
					if (!traversed.contains(linkedBean) && selector.shouldTraverse(linkedBean)) {
						newBeans.add(linkedBean);
					}
				}
			}
			
			for (Dataset newBean : newBeans) {
				traversed.add(newBean);
				conditionallySelect(selector, selected, newBean);

			}
			
		} while (newBeans.size() > 0); // iterate as long as we make progress
		
		return selected;
	}

	/**
	 * @return outgoing links of given type.
	 */
	public List<Dataset> getLinkTargets(Link... types) {
		return getLinkedBeans(types, outgoingLinks);
	}

	/**
	 * @return incoming links of given type.
	 */
	public List<Dataset> getLinkSources(Link... types) {
		return getLinkedBeans(types, incomingLinks);
	}
	
	public URL getContentUrl() {
		return url;
	}


	public void setContentUrl(URL contentUrl) {
		this.url = contentUrl;
	}

	/**
	 * Should only be used during saving or loading.
	 * 
	 * @param handler
	 */
	public DataBeanHandler getHandler() {
		return this.handler;
	}

	/**
	 * Should only be used during saving or loading.
	 * 
	 * @param handler
	 */
	public void setHandler(DataBeanHandler handler) {
		this.handler = handler;
	}



	public DataBeanType getType() {
		return type;
	}



	public void setType(DataBeanType type) {
		this.type = type;
	}



	public String getRepositoryName() {
		if (this.repositoryName == null) {
			return "";
		} else {
			return this.repositoryName;
		}
	}


	public void setRepositoryName(String repositoryName) {
		this.repositoryName = repositoryName;
	}

	/**
	 * Indicate whether the contents have been changed since the contents
	 * were last time uploaded to file broker.
	 * 
	 * FIXME Think about locking.
	 * 
	 */
	public boolean isContentChanged() {
		return this.contentChanged;
	}



	/**
	 * Set content changed status. Should be called with true every time
	 * content is changed.
	 * 
	 * FIXME Think about locking.
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



	public ReentrantReadWriteLock getLock() {
		return this.lock;
	}



	/**
	 * @return A String presentation of this dataset (in this case,
	 * 		   simply the name, to be shown on e.g. the tree).
	 */
	public String toString() {
		return getName();
	}



	private void conditionallySelect(DataBeanSelector selector, LinkedList<Dataset> selected, Dataset bean) {
		if (!selected.contains(bean) && selector.shouldSelect(bean)) {
			selected.add(bean);
		}
	}



	private List<Dataset> getLinkedBeans(Link[] types, LinkedList<LinkedBean> links) {
		LinkedList<Dataset> targets = new LinkedList<Dataset>();
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



	private InputStream getRawContentByteStream() throws IOException {
		if (this.handler == null) {
			throw new IllegalStateException("Handler is null.");
		}
		return handler.getInputStream(this);
	}
}