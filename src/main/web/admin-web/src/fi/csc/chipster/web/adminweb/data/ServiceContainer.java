package fi.csc.chipster.web.adminweb.data;
import java.io.IOException;
import java.io.Serializable;
import java.net.SocketTimeoutException;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.locks.Lock;
import javax.jms.JMSException;
import com.vaadin.data.util.BeanItemContainer;
import com.vaadin.ui.Notification;
import com.vaadin.ui.Notification.Type;

import fi.csc.chipster.web.adminweb.ChipsterConfiguration;
import fi.csc.chipster.web.adminweb.ui.ServicesView;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.AdminAPI;
import fi.csc.microarray.messaging.AdminAPI.AdminAPILIstener;
import fi.csc.microarray.messaging.AdminAPI.NodeStatus;
import fi.csc.microarray.messaging.AdminAPI.NodeStatus.Status;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.Topics;

@SuppressWarnings("serial")
public class ServiceContainer extends BeanItemContainer<ServiceEntry> implements
Serializable {
	
	
	public static final String NAME = "name";
	public static final String COUNT = "count";
	public static final String HOST = "host";
	public static final String STATUS = "status";

	public static final Object[] NATURAL_COL_ORDER  = new String[] {
		NAME,			HOST, 		STATUS };
	public static final String[] COL_HEADERS_ENGLISH = new String[] {
		"Service name", "Host", 	"Status" };
	
	public static final String[] SERVER_NAMES = new String[] { 
		"authenticator", "analyser", "filebroker", "manager" };

	public ServiceContainer() throws InstantiationException,
	IllegalAccessException {
		super(ServiceEntry.class);
	}

	public void update(final ServicesView view) {
		
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
									
									/* Following operation has to lock table component, because addBean() will 
									 * eventually modify its user interface. Keep the lock during the update loop
									 * to avoid showing inconsistent state during the loop.
									 */
									Lock tableLock = view.getTable().getUI().getSession().getLockInstance();
									tableLock.lock();
									try {
										
										removeAllItems();

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

													addBean(service);
												}
											}
										}
										
										//Add a placeholder for each missing server component
										for (String name : ServiceContainer.SERVER_NAMES) {
											
											if (!contains(name)) {
												ServiceEntry entry = new ServiceEntry();
												entry.setName(name);
												entry.setStatus(Status.UNKNOWN);
												addBean(entry);
											}
										}
									} finally {
										tableLock.unlock();
									}
								}
							});

					//Wait for responses									
					api.areAllServicesUp(true);					
					
					endpoint.close();										
					
					view.updateDone();

				} catch (MicroarrayException e) {
					if (e.getCause() != null) {
						//The cause has better message, at least when the broker is not available 
						Throwable cause = e.getCause();
						Notification notification = new Notification(cause.getClass().getSimpleName(), cause.getMessage(), Type.ERROR_MESSAGE); 
						notification.show(view.getApp().getPage());
					}
					e.printStackTrace();
				} catch (JMSException e) {
					e.printStackTrace();
				} catch (InterruptedException e) {
					e.printStackTrace();
				} catch (IOException e) {
					e.printStackTrace();
				} catch (IllegalConfigurationException e) {
					e.printStackTrace();
				} 
			}
		});
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