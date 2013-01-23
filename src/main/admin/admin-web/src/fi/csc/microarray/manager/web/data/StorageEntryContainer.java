package fi.csc.microarray.manager.web.data;

import java.io.Serializable;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.locks.Lock;

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

		ExecutorService execService = Executors.newCachedThreadPool();
		execService.execute(new Runnable() {

			public void run() {
				
				//Simulate some delay
				try {
					Thread.sleep(2000);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

//				try {
//
//					NodeBase nodeSupport = new NodeBase() {
//						public String getName() {
//							return "chipster-admin-web";
//						}
//					};
//
//					ChipsterConfiguration.init();
//					MessagingEndpoint endpoint = new MessagingEndpoint(nodeSupport);
//					AdminAPI api = new AdminAPI(
//							endpoint.createTopic(Topics.Name.ADMIN_TOPIC, AccessMode.READ), new AdminAPILIstener() {
//
//								public void statusUpdated(Map<String, NodeStatus> statuses) {

									/* Following operation has to lock table component, because addBean() will 
									 * eventually modify its user interface. Keep the lock during the update loop
									 * to avoid showing inconsistent state during the loop.
									 */		
				
									//Following will throw nullPointerException if data loading in this thread
									//was faster than UI initialisation in another thread
									Lock tableLock = view.getEntryTable().getUI().getSession().getLock();
									tableLock.lock();
									try {
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

											addBean(entry);
										}

									} finally {
										tableLock.unlock();
									}
//								}
//							});
//
//					//Wait for responses									
//					e.g. api.areAllServicesUp(true);					
//
//					endpoint.close();										

					view.entryUpdateDone();

//				} catch (MicroarrayException e) {
//					e.printStackTrace();
//				} catch (JMSException e) {
//					e.printStackTrace();
//				} catch (InterruptedException e) {
//					e.printStackTrace();
//				} 
			}
		});
	}

	public void showUser(String username) {

		this.removeContainerFilters(USERNAME);

		if (username != null) {
			this.addContainerFilter(USERNAME, username, false, true);
		}
	}
}