package fi.csc.microarray.databeans;

import java.net.URL;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.databeans.DataManager.ContentLocation;
import fi.csc.microarray.databeans.DataManager.StorageMethod;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.QueryResult;
import fi.csc.microarray.databeans.features.RequestExecuter;
import fi.csc.microarray.security.CryptoKey;

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
 * <p>The complete type of a bean is made up from the MIME content type and 
 * the set of type tags. No other information should be used in typing
 * the beans, such as when deciding what operations can be applied to 
 * some particular bean.</p>  
 * 
 * @author Aleksi Kallio
 *
 */
public class DataBean extends DataItemBase {
	

	/**
	 * What to return when none of the data locations are available?
	 */
	public static enum DataNotAvailableHandling {
		
		/**
		 * Throw an exception if data is not available.
		 */
		EXCEPTION_ON_NA,

		/**
		 * Return null if data is not available.
		 */
		NULL_ON_NA,
		
		/**
		 * Return empty data if data is not available.
		 */
		EMPTY_ON_NA,
		
		/**
		 * Return an info string describing the situation and giving details on data availability.
		 */
		INFOTEXT_ON_NA
	}
	
	/**
	 * Traversal specifies the way of traversing links. 
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

	private DataManager dataManager;
	private HashMap<String, Object> contentBoundCache = new HashMap<String, Object>();
	
	private ReentrantReadWriteLock lock = new ReentrantReadWriteLock(true);

	private LinkedList<LinkedBean> outgoingLinks = new LinkedList<LinkedBean>();
	private LinkedList<LinkedBean> incomingLinks = new LinkedList<LinkedBean>();

	private LinkedList<TypeTag> tags = new LinkedList<TypeTag>();
	private boolean tagsSet = false;
	
	private String id;
	
	/**
	 * Timestamp for creation time.
	 */
	private Date date;
	
	private OperationRecord operationRecord;
	private String notes;
	private LinkedHashMap<String, String> toolVersions = new LinkedHashMap<String, String>();
	
	private ContentType contentType;
	private LinkedList<ContentLocation> contentLocations = new LinkedList<ContentLocation>();
	
	private Long size;
	private String checksum;
	private Integer x;
	private Integer y;
	
	
	public DataBean(String name, ContentType contentType, DataManager manager) {
		this(name, contentType, manager, CryptoKey.generateRandom());
	}

	public DataBean(String name, ContentType contentType, DataManager manager, String dataId) {
		this.name = name;
		this.contentType = contentType;
		this.dataManager = manager;
		this.id = dataId;
		this.date = new Date();
		this.toolVersions.put("Chipster", ApplicationConstants.VERSION);
	}

	public String getId() {
		return this.id;
	}
	
	public OperationRecord getOperationRecord() {
		return operationRecord;
	}


	public void setOperationRecord(OperationRecord operationRecord) {
		this.operationRecord = operationRecord;
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


	/**
	 * Only use where really needed.
	 * 
	 * @param newId
	 */
	public void setId(String newId) {
		this.id = newId;
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


	public void delete() {
//		lock.writeLock().lock();
		try {			
			for (ContentLocation contentLocation : contentLocations) {
				contentLocation.getHandler().markDeletable(contentLocation);
			}
			this.contentType = null;			
		} finally {
//			lock.writeLock().unlock();
		}
	}



	/**
	 * Puts named object to content bound cache. Content bound cache is a hashmap type cache where 
	 * objects can be stored. Cache is emptied every time bean content is changed, so it is suited
	 * for caching results that are derived from contents of a single bean. Cache is not persistent,
	 * and generally, user should never assume cached values to be found.  
	 */
	public void putToContentBoundCache(String name, Object value) {
		this.contentBoundCache.put(name, value);
	}



	/**
	 * Gets named object from content cache.
	 * 
	 * @see #addToContentCache(String, Object)
	 * @param name
	 * @return
	 */
	public Object getFromContentBoundCache(String name) {
		return this.contentBoundCache.get(name);
	}


	protected void resetContentBoundCache() {
		this.contentBoundCache.clear();
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
	 * @param storageMethod
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
	
	/**
	 * Attaches a type tag to this data bean, if not already
	 * attached. Type tag represents an aspect of this databean's type.
	 * 
	 * @param tag type tag to be attach
	 * 
	 * @see TypeTag
	 */
	public void addTypeTag(TypeTag tag) {
		if (!tags.contains(tag)) {
			tags.add(tag);
		}
	}
	
	/**
	 * Removes type tag from this databean.
	 * If tag is not attached nothing is done.
	 * 
	 * @param tag tag to be removed
	 * 
	 * @see TypeTag
	 */
	public void removeTypeTag(TypeTag tag) {
		tags.remove(tag);
	}
	
	/**
	 * Gets all type tags attached to this databean.
	 * The set of type tags and MIME content type forms the
	 * complete type of the bean.
	 * 
	 * @return attached type tags
	 * 
	 * @see TypeTag
	 * @see #getContentType()
	 */
	public List<TypeTag> getTypeTags() {
		return tags;
	}
	
	/**
	 * Checks if any of the tags is attached to this data bean.
	 * 
	 * @param alternativeTags type tags to check 
	 * 
	 * @return true iff one or more of the tags are attached
	 * 
	 * @see TypeTag
	 */
	public boolean hasTypeTag(TypeTag... alternativeTags) {
		for (TypeTag alternativeTag : alternativeTags) {
			if (tags.contains(alternativeTag)) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Gives one (arbitrary) location for the content of this DataBean that uses
	 * one of the given storage methods. 
	 * 
	 * @param methods returned ContentLocation must use one of these
	 */
	ContentLocation getContentLocation(StorageMethod... methods) {
		
		if (methods.length == 0) {
			// does not make sense if no methods are specified, check here to avoid programming mistakes
			throw new IllegalArgumentException("must specify at least one method"); 
		}
		
		List<ContentLocation> locations = getContentLocations(methods);
		return locations.isEmpty() ? null : locations.get(0);
	}

	/**
	 * Gives all locations for the content of this DataBean.
	 */
	List<ContentLocation> getContentLocations() {
		
		return contentLocations; 
	}

	/**
	 * Gives all locations for the content of this DataBean that use
	 * one of the given storage methods. 
	 * 
	 * @param methods returned ContentLocations must use one of these
	 */
	List<ContentLocation> getContentLocations(StorageMethod... methods) {
		
		// Collect content locations with matching method
		LinkedList<ContentLocation> locations = new LinkedList<ContentLocation>();
		for (ContentLocation contentLocation : contentLocations) {
			for (StorageMethod method : methods) {
				if (contentLocation.method == method) {
					locations.add(contentLocation);
				}
			}
		}

		return locations; 
	}
	
	/**
	 * Add ContentLocations to bean. Should be used only by DataManager. Others
	 * use DataManager to do this.
	 * 
	 * @see DataManager#addContentLocationForDataBean(DataBean, StorageMethod, URL)
	 * 
	 */
	void addContentLocation(ContentLocation contentLocation) {
		contentLocations.add(contentLocation);
	}
	
	/**
	 * Remove ContentLocation from bean. Should be used only by DataManager. Others
	 * use DataManager to do this.
	 * 
	 * 
	 */
	void removeContentLocation(ContentLocation contentLocation) {
		contentLocations.remove(contentLocation);
	}


	/**
	 * Remove all ContentLocations with given types from bean. Should be used only by DataManager. Others
	 * use DataManager to do this.
	 * 
	 * 
	 */
	void removeContentLocations(StorageMethod... methods) {
		
		// gather
		List<ContentLocation> locationsToBeRemoved = new LinkedList<ContentLocation>();
		for (StorageMethod method : methods) {
			for (ContentLocation location : contentLocations) {
				if (method.equals(location.getMethod())) {
					locationsToBeRemoved.add(location);
				}
			}
		}

		// remove
		for (ContentLocation location : locationsToBeRemoved) {
			contentLocations.remove(location);
		}
	}

	public void setCreationDate(Date date) {
		this.date = date;
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

	/**
	 * Get the oldest unique ancestor the databean was derived from.
	 * 
	 * @param dataBean
	 * @return the given databean if there is no unique ancestor
	 */
	public DataBean getUniqueAncestorRecursively(DataBean dataBean) {
		if (dataBean.getLinkTargets(Link.derivationalTypes()).size() != 1) {
			return dataBean;
		} else {
			return getUniqueAncestorRecursively(dataBean.getLinkTargets(Link.derivationalTypes()).get(0));
		}
	}


	private void conditionallySelect(DataBeanSelector selector, LinkedList<DataBean> selected, DataBean bean) {
		if (!selected.contains(bean) && selector.shouldSelect(bean)) {
			selected.add(bean);
		}
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

	public Long getSize() {
		return size;
	}

	public void setSize(Long size) {
		this.size = size;
	}

	public String getChecksum() {
		return checksum;
	}

	public void setChecksum(String checksum) {
		this.checksum = checksum;
	}

	public boolean isTagsSet() {
		return tagsSet;
	}

	public void setTagsSet(boolean tagsSet) {
		this.tagsSet = tagsSet;
	}

	public void setPosition(Integer x, Integer y) {
		this.x = x;
		this.y = y;
	}

	public Integer getX() {
		return this.x;
	}
	
	public Integer getY() {
		return this.y;
	}

	public LinkedHashMap<String, String> getToolVersions() {
		return this.toolVersions;
	}
}