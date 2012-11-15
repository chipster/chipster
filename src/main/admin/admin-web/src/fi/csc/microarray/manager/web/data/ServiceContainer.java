package fi.csc.microarray.manager.web.data;

import java.io.Serializable;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.jms.JMSException;

import com.vaadin.data.util.BeanItemContainer;

import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.manager.web.ChipsterConfiguration;
import fi.csc.microarray.manager.web.ui.ServicesView;
import fi.csc.microarray.messaging.AdminAPI;
import fi.csc.microarray.messaging.AdminAPI.AdminAPILIstener;
import fi.csc.microarray.messaging.AdminAPI.NodeStatus;
import fi.csc.microarray.messaging.AdminAPI.NodeStatus.Status;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.Topics;

public class ServiceContainer extends BeanItemContainer<ServiceEntry> implements
Serializable {
	
	
	public static final String NAME = "name";
	public static final String COUNT = "count";
	public static final String HOST = "host";
	public static final String STATUS = "status";

	public static final Object[] NATURAL_COL_ORDER  = new String[] {
		NAME, 			COUNT, 				HOST, 		STATUS };
	public static final String[] COL_HEADERS_ENGLISH = new String[] {
		"Service name", "Service count", 	"Host", 	"Status" };
	
	public static final String[] SERVER_NAMES = new String[] { 
		"authenticator", "analyser", "filebroker", "manager" };

	public ServiceContainer() throws InstantiationException,
	IllegalAccessException {
		super(ServiceEntry.class);
	}

	public void update(final ServicesView view) {
		
		//Add a placeholder for each required server component
		for (String name : ServiceContainer.SERVER_NAMES) {
			
			if (!contains(name)) {
				ServiceEntry entry = new ServiceEntry();
				entry.setName(name);
				entry.setStatus(Status.UNKNOWN);
				addBean(entry);
			}
		}
		
		ExecutorService execService = Executors.newCachedThreadPool();
		execService.execute(new Runnable() {

			public void run() {

				try {

					NodeBase nodeSupport = new NodeBase() {
						public String getName() {
							return "chipster-admin-web";
						}
					};

					ChipsterConfiguration.init();
					MessagingEndpoint endpoint = new MessagingEndpoint(nodeSupport);
					AdminAPI api = new AdminAPI(
							endpoint.createTopic(Topics.Name.ADMIN_TOPIC, AccessMode.READ), new AdminAPILIstener() {

								public void statusUpdated(Map<String, NodeStatus> statuses) {

									//removeAllItems();

									for (Entry<String, NodeStatus> entry : statuses.entrySet()) {
										NodeStatus node = entry.getValue();
										
										if (node.host != null) {
											String hosts[] = node.host.split(", ");
											for (String host : hosts) {

												ServiceEntry service = new ServiceEntry();
												service.setName(node.name);
												service.setHost(host);
												service.setStatus(node.status);
												service.setCount(node.count);

												synchronized (view.getApp()) {
													
													removeAll(service.getName());
													
													addBean(service);
													view.dataUpdated();
												}
											}
										}
									}					
								}
							});

					//Wait for responses
					api.areAllServicesUp(true);

					endpoint.close();
			

				} catch (MicroarrayException e) {
					e.printStackTrace();
				} catch (JMSException e) {
					e.printStackTrace();
				} catch (InterruptedException e) {
					e.printStackTrace();
				} 
			}
		});
	}
	
	private void removeAll(String name) {
		
		//Concurrent deletion and iteration not supported 
		List<ServiceEntry> removeList = new LinkedList<ServiceEntry>();
		
		for (ServiceEntry entry : getItemIds()) {
			if (name.equals(entry.getName())) {
				removeList.add(entry);
			}
		}
		
		for (ServiceEntry entry : removeList) {
			removeItem(entry);
		}
	}
	
	private boolean contains(String name) {
		for (ServiceEntry entry : getItemIds()) {
			if (name.equals(entry.getName())) {
				return true;
			}
		}
		return false;
	}
}