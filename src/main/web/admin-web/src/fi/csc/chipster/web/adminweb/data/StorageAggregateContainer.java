package fi.csc.chipster.web.adminweb.data;

import java.io.Serializable;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.locks.Lock;

import javax.jms.JMSException;

import com.vaadin.data.util.BeanItemContainer;

import fi.csc.chipster.web.adminweb.ChipsterConfiguration;
import fi.csc.chipster.web.adminweb.ui.StorageView;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.TempTopicMessagingListenerBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;

public class StorageAggregateContainer extends BeanItemContainer<StorageAggregate> implements
Serializable {
	
	public static final String USERNAME = "username";
	public static final String SIZE = "size";

	public static final Object[] NATURAL_COL_ORDER  = new String[] {
		USERNAME, 		SIZE };

	public static final String[] COL_HEADERS_ENGLISH = new String[] {
		"Username", 	"Total size" };
	
	

	public final String TOTAL_USERNAME = "TOTAL";

//	private StorageEntryContainer entryContainer;

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
		
//		this.entryContainer = entryContainer;
	}

//	public StorageAggregate update(final StorageView view) {
//
//		removeAllItems();
//		
//		long totalSize = 0;
//
//		HashMap<String, Long> aggregateMap = new HashMap<String, Long>();
//		
//		entryContainer.removeAllContainerFilters();
//		
//		for (StorageEntry entry : entryContainer.getItemIds()) {
//			
//			String username = entry.getUsername();
//			if (!aggregateMap.containsKey(username)) {
//				aggregateMap.put(username, entry.getSize());
//			} else {
//				aggregateMap.put(username, aggregateMap.get(username) + entry.getSize());
//			}
//			
//			totalSize += entry.getSize();
//		}
//		
//		StorageAggregate totalBean = new StorageAggregate();
//		totalBean.setUsername(TOTAL_USERNAME);
//		totalBean.setSize(totalSize);
//		this.addBean(totalBean);
//		
//		for (Entry<String, Long> aggregate : aggregateMap.entrySet()) {
//			
//			StorageAggregate bean = new StorageAggregate();
//			bean.setUsername(aggregate.getKey());
//			bean.setSize(aggregate.getValue());
//			
//			this.addBean(bean);
//		}
//				
//		this.diskUsage = totalSize;
//		this.diskAvailable = 500000000000l;
//		
//		return totalBean;
//	}
	
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

					CommandMessage request = new CommandMessage(CommandMessage.COMMAND_LIST_STORAGE_USAGE_OF_USERS);
					final CountDownLatch latch = new CountDownLatch(1);

					// TODO clean up this
					StorageAggregateMessageListener replyListener = new StorageAggregateMessageListener(latch);
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

							for (StorageAggregate entry : replyListener.getEntries()) {
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
	
	
	private class StorageAggregateMessageListener extends TempTopicMessagingListenerBase {
		
		private CountDownLatch latch;
		private List<StorageAggregate> entries;
		
		public StorageAggregateMessageListener(CountDownLatch latch) {
			this.latch = latch;
		}
		
		public void onChipsterMessage(ChipsterMessage msg) {
			ParameterMessage resultMessage = (ParameterMessage) msg;
			
			String namesString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_USERNAME_LIST);
			String sizesString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_SIZE_LIST);
			
			String[] names = namesString.split("\t");
			String[] sizes = sizesString.split("\t");
			
			entries = new LinkedList<StorageAggregate>();
			for (int i = 0; i < names.length && i < sizes.length; i++) {
				
				StorageAggregate entry = new StorageAggregate();
				entry.setUsername(names[i]);
				entry.setSize(Long.parseLong(sizes[i]));				
				entries.add(entry);
			}
								
			latch.countDown();
		}

		public List<StorageAggregate> getEntries() {
			return entries;
		}
		
	}
}