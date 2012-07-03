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

	private StorageEntryContainer enryContainer;

	public StorageAggregateContainer(StorageEntryContainer entryContainer) throws InstantiationException,
	IllegalAccessException {
		super(StorageAggregate.class);
		
		this.enryContainer = entryContainer;
	}

	public StorageAggregate update(final StorageView view) {

		removeAllItems();
		
		long totalSize = 0;

		HashMap<String, Long> aggregateMap = new HashMap<String, Long>();

		for (StorageEntry entry : enryContainer.getItemIds()) {
			
			String username = entry.getUsername();
			if (!aggregateMap.containsKey(username)) {
				aggregateMap.put(username, entry.getSize());
			} else {
				aggregateMap.put(username, aggregateMap.get(username) + entry.getSize());
			}
		}
		
		for (Entry<String, Long> aggregate : aggregateMap.entrySet()) {
			
			StorageAggregate bean = new StorageAggregate();
			bean.setUsername(aggregate.getKey());
			bean.setSize(aggregate.getValue());
			
			this.addBean(bean);
			
			totalSize += bean.getSize();
		}
		
		StorageAggregate totalBean = new StorageAggregate();
		totalBean.setUsername(TOTAL_USERNAME);
		totalBean.setSize(totalSize);
		
		this.addBean(totalBean);
		
		return totalBean;
	}
}