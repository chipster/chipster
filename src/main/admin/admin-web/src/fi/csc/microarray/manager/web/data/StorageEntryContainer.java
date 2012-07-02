package fi.csc.microarray.manager.web.data;

import java.io.Serializable;
import java.util.Map;
import java.util.Map.Entry;

import javax.jms.JMSException;

import com.vaadin.data.util.BeanItemContainer;

import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.manager.web.ui.ServicesView;
import fi.csc.microarray.messaging.AdminAPI;
import fi.csc.microarray.messaging.AdminAPI.AdminAPILIstener;
import fi.csc.microarray.messaging.AdminAPI.NodeStatus;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.Topics;

public class StorageEntryContainer extends BeanItemContainer<Service> implements
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
		super(Service.class);
	}

	public void update(final ServicesView view) {

		new Runnable() {

			public void run() {

				try {

					NodeBase nodeSupport = new NodeBase() {
						public String getName() {
							return "chipster-admin-web";
						}
					};

					MessagingEndpoint endpoint = new MessagingEndpoint(nodeSupport);
					AdminAPI api = new AdminAPI(
							endpoint.createTopic(Topics.Name.ADMIN_TOPIC, AccessMode.READ), new AdminAPILIstener() {

								public void statusUpdated(Map<String, NodeStatus> statuses) {

									removeAllItems();


									for (Entry<String, NodeStatus> entry : statuses.entrySet()) {
										NodeStatus node = entry.getValue();
										
										String hosts[] = node.host.split(", ");
										for (String host : hosts) {

											Service service = new Service();
											service.setName(node.name);
											service.setHost(host);
											service.setStatus(node.status);
											service.setCount(node.count);

											addBean(service);
										}
									}					
								}
							});

					//Wait for responses
					api.areAllServicesUp(true);

					endpoint.close();
					
					synchronized (view.getApp()) {
						view.dataUpdated();
					}

				} catch (MicroarrayException e) {
					e.printStackTrace();
				} catch (JMSException e) {
					e.printStackTrace();
				} catch (InterruptedException e) {
					e.printStackTrace();
				} 
			}
		}.run();
	}
}