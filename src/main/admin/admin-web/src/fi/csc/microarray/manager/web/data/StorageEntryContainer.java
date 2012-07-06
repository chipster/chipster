package fi.csc.microarray.manager.web.data;

import java.io.Serializable;
import java.util.Random;

import com.vaadin.data.util.BeanItemContainer;

import fi.csc.microarray.manager.web.ui.StorageView;
import fi.csc.microarray.manager.web.util.RandomUtil;

public class StorageEntryContainer extends BeanItemContainer<StorageEntry> implements
Serializable {

	/**
	 * Natural property order for Service bean. Used in tables and forms.
	 */
	public static final Object[] NATURAL_COL_ORDER = new Object[] {
		"username", "name", "size", "date", "deleteLink" };

	/**
	 * "Human readable" captions for properties in same order as in
	 * NATURAL_COL_ORDER.
	 */
	public static final String[] COL_HEADERS_ENGLISH = new String[] {
		"Username", "Session name", "Size", "Date", " "};

	public StorageEntryContainer() throws InstantiationException,
	IllegalAccessException {
		super(StorageEntry.class);
	}

	public void update(final StorageView view) {
		
		final int COUNT = 300;
		
		removeAllItems();
		
    	
        Random rnd = new Random();

        StorageEntry entry;

        for (int i = 0; i < COUNT; i++) {
            entry = new StorageEntry();
            
            entry.setDate(RandomUtil.getRandomDate(rnd, 2011));
            entry.setUsername(RandomUtil.getRandomUserName(rnd));
            entry.setSize(Math.abs(rnd.nextInt(rnd.nextInt(9000000))*1000l));
            entry.setName(RandomUtil.getRandomSessionName(rnd));
            
            this.addBean(entry);
        }
    }

	public void showUser(String username) {
		
		this.removeContainerFilters("username");
		
		if (username != null) {
			this.addContainerFilter("username", username, false, true);
		}
	}
}