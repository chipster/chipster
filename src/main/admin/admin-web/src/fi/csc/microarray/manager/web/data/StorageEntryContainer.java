package fi.csc.microarray.manager.web.data;

import java.io.Serializable;
import java.util.Date;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.locks.Lock;

import javax.jms.JMSException;

import com.vaadin.data.util.BeanItemContainer;

import fi.csc.microarray.manager.web.ChipsterConfiguration;
import fi.csc.microarray.manager.web.ui.StorageView;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.TempTopicMessagingListenerBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;

@SuppressWarnings("serial")
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

				
				MessagingEndpoint endpoint = null;
				try {

					NodeBase nodeSupport = new NodeBase() {
						public String getName() {
							return "admin";
						}
					};

					ChipsterConfiguration.init();
					endpoint = new MessagingEndpoint(nodeSupport);
					
					// TODO close topic
					MessagingTopic filebrokerAdminTopic = endpoint.createTopic(Topics.Name.FILEBROKER_ADMIN_TOPIC, AccessMode.WRITE);

					CommandMessage request = new CommandMessage(CommandMessage.COMMAND_GET_SESSIONS_FOR_USER);
					request.addNamedParameter("username", "testiuuseri");
					final CountDownLatch latch = new CountDownLatch(1);

					// TODO clean up this
					StorageEntryMessageListener replyListener = new StorageEntryMessageListener(latch);
					filebrokerAdminTopic.sendReplyableMessage(request, replyListener);

					// wait for responses TODO timeout									
					latch.await();

					// TODO check if results, timeout

					/* Following operation has to lock table component, because addBean() will 
					 * eventually modify its user interface. Keep the lock during the update loop
					 * to avoid showing inconsistent state during the loop.
					 */		

					//Following is null if data loading in this thread
					//was faster than UI initialisation in another thread
if (view.getEntryTable().getUI() != null) {
	Lock tableLock = view.getEntryTable().getUI().getSession().getLockInstance();
					tableLock.lock();
					try {
						removeAllItems();

						for (StorageEntry entry : replyListener.getEntries()) {
							addBean(entry);
						}

					} finally {
						tableLock.unlock();
					}
}
					
					// TODO should be in the last finally?
					view.entryUpdateDone();

				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} finally {
					if (endpoint != null) {
						try {
							endpoint.close();
						} catch (JMSException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
					
				}
			}});
	}

	public void showUser(String username) {

		this.removeContainerFilters(USERNAME);

		if (username != null) {
			this.addContainerFilter(USERNAME, username, false, true);
		}
	}

	

	
	
	private class StorageEntryMessageListener extends TempTopicMessagingListenerBase {
		
		private CountDownLatch latch;
		private List<StorageEntry> entries;
		
		public StorageEntryMessageListener(CountDownLatch latch) {
			this.latch = latch;
		}
		
		public void onChipsterMessage(ChipsterMessage msg) {
			ParameterMessage resultMessage = (ParameterMessage) msg;
			
			entries = new LinkedList<StorageEntry>();
			StorageEntry entry = new StorageEntry();
			entry.setDate(new Date());
			entry.setUsername(resultMessage.getNamedParameter("neppi"));
			entry.setSize(122);
			entry.setName("sessssio");
			entries.add(entry);
			
			latch.countDown();
		}

		public List<StorageEntry> getEntries() {
			return entries;
		}
		
	}



}