package fi.csc.microarray.manager.web.data;

import java.io.Serializable;
import java.util.Random;

import com.vaadin.data.hbnutil.ContainerFilter;
import com.vaadin.data.hbnutil.StringContainerFilter;
import com.vaadin.data.util.BeanItemContainer;

import fi.csc.microarray.manager.web.ui.StorageView;

public class StorageEntryContainer extends BeanItemContainer<StorageEntry> implements
Serializable {

	/**
	 * Natural property order for Service bean. Used in tables and forms.
	 */
	public static final Object[] NATURAL_COL_ORDER = new Object[] {
		"username", "name", "size", "date" };

	/**
	 * "Human readable" captions for properties in same order as in
	 * NATURAL_COL_ORDER.
	 */
	public static final String[] COL_HEADERS_ENGLISH = new String[] {
		"Username", "Session name", "Size", "Date" };

	public StorageEntryContainer() throws InstantiationException,
	IllegalAccessException {
		super(StorageEntry.class);
	}

	public void update(final StorageView view) {
		
		final int COUNT = 100;
		
		removeAllItems();
		
    	String[] entryName = new String[] { "ngs-session1", "ngs-session2", "ngs-session3" };
    	
        Random rnd = new Random();

        StorageEntry entry;

        for (int i = 0; i < COUNT; i++) {
            entry = new StorageEntry();
            
            entry.setDate(RandomUtil.getRandomDate(rnd));
            entry.setUsername(RandomUtil.getRandomUserName(rnd));
            entry.setSize(Math.abs(rnd.nextInt(20000000)*1000l));
            entry.setName(entryName[rnd.nextInt(entryName.length)]);
            
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