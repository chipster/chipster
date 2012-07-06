package fi.csc.microarray.manager.web.data;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map.Entry;

import com.vaadin.data.util.BeanItemContainer;

import fi.csc.microarray.manager.web.ui.StorageView;

public class StorageAggregateContainer extends BeanItemContainer<StorageAggregate> implements
Serializable {

	/**
	 * Natural property order for Service bean. Used in tables and forms.
	 */
	public static final Object[] NATURAL_COL_ORDER = new Object[] {
		"username", "size" };

	/**
	 * "Human readable" captions for properties in same order as in
	 * NATURAL_COL_ORDER.
	 */
	public static final String[] COL_HEADERS_ENGLISH = new String[] {
		"Username", "Total size" };

	public final String TOTAL_USERNAME = "TOTAL";

	private StorageEntryContainer entryContainer;

	public long getDiskUsage() {
		return diskUsage;
	}

	public long getDiskAvailable() {
		return diskAvailable;
	}

	private long diskUsage = 0;
	private long diskAvailable = 0;

	public StorageAggregateContainer(StorageEntryContainer entryContainer) throws InstantiationException,
	IllegalAccessException {
		super(StorageAggregate.class);
		
		this.entryContainer = entryContainer;
	}

	public StorageAggregate update(final StorageView view) {

		removeAllItems();
		
		long totalSize = 0;

		HashMap<String, Long> aggregateMap = new HashMap<String, Long>();
		
		entryContainer.removeAllContainerFilters();
		
		for (StorageEntry entry : entryContainer.getItemIds()) {
			
			String username = entry.getUsername();
			if (!aggregateMap.containsKey(username)) {
				aggregateMap.put(username, entry.getSize());
			} else {
				aggregateMap.put(username, aggregateMap.get(username) + entry.getSize());
			}
			
			totalSize += entry.getSize();
		}
		
		StorageAggregate totalBean = new StorageAggregate();
		totalBean.setUsername(TOTAL_USERNAME);
		totalBean.setSize(totalSize);
		this.addBean(totalBean);
		
		for (Entry<String, Long> aggregate : aggregateMap.entrySet()) {
			
			StorageAggregate bean = new StorageAggregate();
			bean.setUsername(aggregate.getKey());
			bean.setSize(aggregate.getValue());
			
			this.addBean(bean);
		}
				
		this.diskUsage = totalSize;
		this.diskAvailable = 500000000000l;
		
		return totalBean;
	}
}