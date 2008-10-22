package fi.csc.microarray.databeans;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Date;
import java.util.List;

import javax.jms.JMSException;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.databeans.features.QueryResult;
import fi.csc.microarray.messaging.message.PayloadMessage;

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
public interface DataBean extends DataItem {
	
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
	
	/**
	 * Returns result of the feature query. Acts as a gateway to Feature API which
	 * is the way bean contents can be used besides using the raw binary content.
	 * Available features depend on what feature factories are plugged to DataManager.
	 */
	public QueryResult queryFeatures(String featureName);
	
	/**
	 * @return The date and time when this dataset was created.
	 */
	public Date getDate();

	/**
	 * @see #getProperty(String)
	 */
	public void setNotes(String note);
	
	public String getNotes();
	
	/**
	 * Creates a link between this and target bean. Links represent relationships
	 * between beans. Links have hardcoded types.
	 * 
	 * @see Link
	 */
	public void addLink(Link type, DataBean target);
	
	/**
	 * Removes a link. Does not remove links of different type between same beans.
	 */
	public void removeLink(Link type, DataBean target);
	
	/**
	 * @return outgoing links of given type.
	 */
	public List<DataBean> getLinkTargets(Link... types);
	
	/**
	 * @return incoming links of given type.
	 */
	public List<DataBean> getLinkSources(Link... types);
	
	
	/**
	 * Convenience method which includes all beans into result.
	 * 
	 * @see #traverseLinks(fi.csc.microarray.databeans.DataBean.Link, fi.csc.microarray.databeans.DataBean.Traversal, DataBeanSelector)
	 */
	public List<DataBean> traverseLinks(Link[] types, Traversal traversal);
	
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
	public List<DataBean> traverseLinks(Link[] types, Traversal traversal, DataBeanSelector selector);
	
	/**
	 * @return The operation that has been selected for this dataset (and which
	 * 		   may have already been done at least once to produce another
	 * 		   dataset, or which has not yet been conducted, and maybe never
	 * 		   will, depending on the user).
	 */
	public Operation getOperation();

	/**
	 * Associates the given operation with this DataBean.
	 * 
	 * @param operation to associate
	 */
	public void setOperation(Operation operation);

	/**
	 * A convenience method for gathering streamed binary content into
	 * a byte array.
	 * 
	 *   @see #getContentByteStream()
	 */
	public byte[] getContents() throws MicroarrayException;

	/**
	 * Returns raw binary content of this bean. Use Feature API for 
	 * higher level accessing.
	 * 
	 * @see #queryFeatures(String)
	 */
	public InputStream getContentByteStream() throws MicroarrayException, IOException;

	/**
	 * Allows rewriting of raw bean content. Close the stream by calling closeContentOutputStream(...)
	 * on the same bean.
	 * @see #closeContentOutputStreamAndUnlockDataBean(OutputStream)
	 */
	public OutputStream getContentOutputStreamAndLockDataBean() throws MicroarrayException, IOException;

	/**
	 * Closes output stream and generates required events.
	 */
	public void closeContentOutputStreamAndUnlockDataBean(OutputStream out) throws MicroarrayException, IOException;

	/**
	 * @return MIME content type of the bean.
	 */
	public ContentType getContentType();
	
	/**
	 * @see #getContentType()
	 */
	public void setContentType(ContentType contentType);
	
	/**
	 * Checks if bean's content type is any of the given alternatives. Case insensitive.
	 */
	public boolean isContentTypeCompatitible(String... alternatives);
	
	// TODO should be integrated and hidden away
	public void initialiseStreamStartCache() throws IOException, MicroarrayException;
	
	public void updateRemoteCache(String payloadName, PayloadMessage payloadMessage) throws JMSException, MicroarrayException, IOException;
	
	/**
	 * Returns content size in bytes.
	 */
	public long getContentLength();
	
	
	/**
	 * Puts named object to content cache. Content cache is a hashmap type per DataBean cache where 
	 * objects can be stored. Cache is emptied every time bean content is changed, so it is suited
	 * for caching results that are derived from contents of a single bean. Cache is not persistant,
	 * and generally, user should never assume cached values to be found.  
	 */
	public void putToContentCache(String name, Object value);

	/**
	 * Gets named object from content cache.
	 * 
	 * @see #addToContentCache(String, Object)
	 * @param name
	 * @return
	 */
	public Object getFromContentCache(String name);
}