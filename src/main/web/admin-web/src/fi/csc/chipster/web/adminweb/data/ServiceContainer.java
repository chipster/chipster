package fi.csc.chipster.web.adminweb.data;
import java.io.Serializable;
import java.util.Map;
import java.util.Map.Entry;

import javax.jms.JMSException;

import com.vaadin.data.util.BeanItemContainer;

import fi.csc.chipster.web.adminweb.ui.ServicesView;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.admin.AdminAPI;
import fi.csc.microarray.messaging.admin.AdminAPI.AdminAPILIstener;
import fi.csc.microarray.messaging.admin.AdminAPI.NodeStatus;
import fi.csc.microarray.messaging.admin.AdminAPI.NodeStatus.Status;

@SuppressWarnings("serial")
public class ServiceContainer extends BeanItemContainer<ServiceEntry> implements AdminAPILIstener, Serializable {
	
	
	public static final String NAME = "name";
	public static final String COUNT = "count";
	public static final String HOST = "host";
	public static final String STATUS = "status";

	public static final Object[] NATURAL_COL_ORDER  = new String[] {
		NAME,			HOST, 		STATUS };
	public static final String[] COL_HEADERS_ENGLISH = new String[] {
		"Service name", "Host", 	"Status" };
	
	public static final String[] SERVER_NAMES = new String[] { 
		"authenticator", "analyser", "filebroker", "manager", "jobmanager" };
	private ServicesView view;

	public ServiceContainer(ServicesView view) throws InstantiationException, IllegalAccessException {
		super(ServiceEntry.class);
		this.view = view;
	}

	public void update() {
		
		try {

			MessagingEndpoint endpoint = view.getApp().getEndpoint();
			AdminAPI api = new AdminAPI(
					endpoint.createTopic(Topics.Name.ADMIN_TOPIC, AccessMode.READ), this);
			
			//Wait for responses									
			api.areAllServicesUp(true);														

		} catch (JMSException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		} 	
	}
	
	public void statusUpdated(Map<String, NodeStatus> statuses) {
		updateUI(view, statuses);
	}

	
	private void updateUI(final ServicesView view, final Map<String, NodeStatus> statuses) {
		/* Following operation has to lock table component, because addBean() will 
		 * eventually modify its user interface. Keep the lock during the update loop
		 * to avoid showing inconsistent state during the loop.
		 */
		view.updateUI(new Runnable() {
			@Override
			public void run() {

				removeAllItems();

				for (Entry<String, NodeStatus> entry : statuses.entrySet()) {

					NodeStatus node = entry.getValue();

					for (String host : node.getHosts()) {

						ServiceEntry service = new ServiceEntry();
						service.setName(node.name);
						service.setHost(host);
						service.setStatus(node.status);
						service.setCount(node.getCount());

						addBean(service);
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