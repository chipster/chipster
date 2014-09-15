package fi.csc.chipster.web.adminweb.data;

import java.io.Serializable;
import java.util.List;
import java.util.concurrent.locks.Lock;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import com.vaadin.data.util.BeanItemContainer;
import com.vaadin.ui.Notification;
import com.vaadin.ui.Notification.Type;

import fi.csc.chipster.web.adminweb.ui.StorageView;
import fi.csc.microarray.messaging.admin.StorageAdminAPI;
import fi.csc.microarray.messaging.admin.StorageAggregate;

public class StorageAggregateContainer extends BeanItemContainer<StorageAggregate> implements Serializable {
	
	private static final Logger logger = Logger.getLogger(StorageAggregateContainer.class);
	
	public static final String USERNAME = "username";
	public static final String SIZE = "size";

	public static final Object[] NATURAL_COL_ORDER  = new String[] {
		USERNAME, 		SIZE };

	public static final String[] COL_HEADERS_ENGLISH = new String[] {
		"Username", 	"Total size" };

	public long getDiskUsage() {
		return diskUsage;
	}

	public long getDiskAvailable() {
		return diskAvailable;
	}

	private long diskUsage = 0;
	private long diskAvailable = 0;
	private StorageAdminAPI adminEndpoint;

	public StorageAggregateContainer(StorageAdminAPI adminEndpoint) throws InstantiationException,
	IllegalAccessException {
		super(StorageAggregate.class);
		this.adminEndpoint = adminEndpoint;
	}
	
	public void update(final StorageView view) {
		
		List<StorageAggregate> entries;
		try {
			entries = adminEndpoint.listStorageUsageOfUsers();		

			if (entries != null) {
				//Following is null if data loading in this thread
				//was faster than UI initialisation in another thread
				if (view.getEntryTable().getUI() != null) {
					Lock tableLock = view.getEntryTable().getUI().getSession().getLockInstance();
					tableLock.lock();
					try {
						removeAllItems();

						for (StorageAggregate entry : entries) {
							addBean(entry);
						}

					} finally {
						tableLock.unlock();
					}
				}		
			} else {
				Notification.show("Timeout", "Chipster filebroker server doesn't respond", Type.ERROR_MESSAGE);
				logger.error("timeout while waiting storage usage of users");
			}
			
		} catch (JMSException | InterruptedException e) {
			logger.error(e);
		}			
	}
}