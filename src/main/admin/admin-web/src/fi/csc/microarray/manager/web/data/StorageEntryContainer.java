package fi.csc.microarray.manager.web.data;

import java.io.Serializable;
import java.util.Random;

import com.vaadin.data.util.BeanItemContainer;

import fi.csc.microarray.manager.web.ui.StorageView;
import fi.csc.microarray.manager.web.util.RandomUtil;

public class StorageEntryContainer extends BeanItemContainer<StorageEntry> implements
Serializable {

	public static final String USERNAME = "username";
	public static final String NAME = "name";
	public static final String SIZE = "size";
	public static final String DATE = "date";
	public static final String DELETE_LINK = "deleteLink";


	public static final Object[] NATURAL_COL_ORDER  = new String[] {
		USERNAME, 		NAME, 			SIZE, 	DATE, 		DELETE_LINK };

	public static final String[] COL_HEADERS_ENGLISH = new String[] {
		"Username", 	"Session name", "Size", "Date", 	" " };

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
		
		this.removeContainerFilters(USERNAME);
		
		if (username != null) {
			this.addContainerFilter(USERNAME, username, false, true);
		}
	}
}