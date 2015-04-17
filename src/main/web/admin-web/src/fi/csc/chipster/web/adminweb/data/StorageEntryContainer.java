package fi.csc.chipster.web.adminweb.data;

import java.io.Serializable;
import java.util.LinkedList;
import java.util.List;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import com.vaadin.data.util.BeanItemContainer;
import com.vaadin.ui.Notification;
import com.vaadin.ui.Notification.Type;

import fi.csc.chipster.web.adminweb.ui.StorageView;
import fi.csc.microarray.messaging.admin.StorageAdminAPI;
import fi.csc.microarray.messaging.admin.StorageEntry;

@SuppressWarnings("serial")
public class StorageEntryContainer extends BeanItemContainer<StorageEntry> implements Serializable {
	
	private static final Logger logger = Logger.getLogger(StorageEntryContainer.class);

	public static final String USERNAME = "username";
	public static final String NAME = "name";
	public static final String SIZE = "size";
	public static final String DATE = "date";
	public static final String DELETE_LINK = "deleteLink";


	public static final Object[] NATURAL_COL_ORDER  = new String[] {
		USERNAME, 		NAME, 			SIZE, 	DATE, 		DELETE_LINK };

	public static final String[] COL_HEADERS_ENGLISH = new String[] {
		"Username", 	"Session name", "Size", "Last access date", 	" " };
	
	private StorageAdminAPI adminEndpoint;
	
	public StorageEntryContainer(StorageAdminAPI adminEndpoint) throws InstantiationException, IllegalAccessException {
		super(StorageEntry.class);
		this.adminEndpoint = adminEndpoint;
	}

	public void update(final StorageView view, final String username) {		
			
		List<StorageEntry> entries;
		try {
			if (username != null) {
				entries = adminEndpoint.listStorageUsageOfSessions(username);
			} else {
				//clear table
				entries = new LinkedList<>();
			}
			
			if (entries != null) {

				updateUI(view, entries);
				
			} else {
				Notification.show("Timeout", "Chipster filebroker server doesn't respond", Type.ERROR_MESSAGE);
				logger.error("timeout while waiting storage usage of sessions");
			}
		} catch (JMSException | InterruptedException e) {
			logger.error(e);
		}
	}

	private void updateUI(StorageView view, final List<StorageEntry> entries) {
		view.updateUI(new Runnable() {
			@Override
			public void run() {
				removeAllItems();

				for (StorageEntry entry : entries) {
					addBean(entry);
				}
			}
		});
	}
}
