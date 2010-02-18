package fi.csc.microarray.databeans;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import org.apache.log4j.Logger;

import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.QueryResult;
import fi.csc.microarray.databeans.features.RequestExecuter;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.InputStreamSource;
import fi.csc.microarray.util.StreamStartCache;
import fi.csc.microarray.util.Strings;

public abstract class DataBeanBase implements DataBean {
	
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
	
	public String toStringRecursively(int i) {
		return Strings.repeat("  ", i) + toString();
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
	
	
}
