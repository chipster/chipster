package fi.csc.chipster.web.adminweb;

import java.io.IOException;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import com.vaadin.annotations.Push;
import com.vaadin.shared.communication.PushMode;
import com.vaadin.shared.ui.ui.Transport;
import com.vaadin.annotations.Theme;
import com.vaadin.server.ClientConnector.DetachListener;
import com.vaadin.server.VaadinRequest;
import com.vaadin.ui.Component;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.Notification;
import com.vaadin.ui.Notification.Type;
import com.vaadin.ui.UI;
import com.vaadin.ui.VerticalLayout;

import fi.csc.chipster.web.adminweb.ui.JobLogView;
import fi.csc.chipster.web.adminweb.ui.JobsView;
import fi.csc.chipster.web.adminweb.ui.ReportView;
import fi.csc.chipster.web.adminweb.ui.ServicesView;
import fi.csc.chipster.web.adminweb.ui.StatView;
import fi.csc.chipster.web.adminweb.ui.StorageView;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.JMSMessagingEndpoint;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.admin.ManagerConfiguration;

@SuppressWarnings("serial")
@Theme("admin")
/*
 * Workaround for Vaadin Push and session-per-request -pattern used by Hibernate
 * HbnContainer's ServletFilter doesn't work with WebSockets
 * (PushMode.AUTOMATIC)
 */ 
@Push(value=PushMode.MANUAL, transport=Transport.LONG_POLLING)
public class ChipsterAdminUI extends UI implements DetachListener {
	
	private static final Logger logger = Logger.getLogger(ChipsterAdminUI.class);

	private HorizontalLayout horizontalSplit;

	private NavigationMenu navigationLayout;;

	private ServicesView serviceView;
	private StorageView storageView;
	private JobsView jobsView;
	private JobLogView jobLogView;
	private StatView statView;
	private ReportView reportView;

	private VerticalLayout emptyView;

	private HorizontalLayout toolbarLayout;

	private JMSMessagingEndpoint endpoint;

	private void buildMainLayout() {

		this.getPage().setTitle("Chipster admin");
		horizontalSplit = new HorizontalLayout();
		horizontalSplit.setSizeFull();

		setContent(horizontalSplit);
		horizontalSplit.setHeight(100, Unit.PERCENTAGE);
		showEmtpyView();	
	}

	private Component getNavigation() {

		if (navigationLayout == null) {
			navigationLayout = new NavigationMenu(this);
			navigationLayout.showDefaultView();
		}
		return navigationLayout;
	}

	private void setMainComponent(Component c) {

		horizontalSplit.removeAllComponents();

		horizontalSplit.addComponent(getNavigation(), 0);
		horizontalSplit.addComponent(c, 1);
		
		horizontalSplit.setExpandRatio(c, 1);
	}

	protected void showServicesView() {
		if (serviceView == null) {
			serviceView = new ServicesView(this);
		}
		setMainComponent(serviceView);
		serviceView.update();
	}

	public void showJobLogView() {
		if (jobLogView == null) {
			jobLogView = new JobLogView(this);
		}
		setMainComponent(jobLogView);
		jobLogView.update();
	}


	public void showJobsView() {
		if (jobsView == null) {
			jobsView = new JobsView(this);
		}
		setMainComponent(jobsView);
		jobsView.update();
	}


	public void showStorageView() {
		if (storageView == null) {
			storageView = new StorageView(this);
		}
		setMainComponent(storageView);
		storageView.update();
	}
	
	public void showStatView() {
		if (statView == null) {
			statView = new StatView(this);
		}
		setMainComponent(statView);
		statView.update();
	}
	
	public void showReportView() {
		if (reportView == null) {
			reportView = new ReportView(this);
		}
		setMainComponent(reportView);
		reportView.update();
	}

	public void showEmtpyView() {
		if (emptyView == null) {
			emptyView = new VerticalLayout();
			
			emptyView.addComponent(getToolbar());
			
			emptyView.setSizeFull();
			emptyView.setStyleName("empty-view");
		}
		setMainComponent(emptyView);
	}
	
	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			
			toolbarLayout = new HorizontalLayout();
			
			Label spaceEater = new Label(" ");
			toolbarLayout.addComponent(spaceEater);
			toolbarLayout.setExpandRatio(spaceEater, 1);
			
			toolbarLayout.addComponent(getTitle());	
			
			toolbarLayout.setWidth("100%");
			toolbarLayout.setHeight(40, Unit.PIXELS);
			toolbarLayout.setStyleName("toolbar");
		}

		return toolbarLayout;
	}
	
	public Component getTitle() {
		Label label = new Label("Chipster admin");
		label.addStyleName("title");
		label.setWidth(250, Unit.PIXELS);
		return label;
	}

	@Override
	protected void init(VaadinRequest request) {
		try {

			NodeBase nodeSupport = new NodeBase() {
				public String getName() {
					return "chipster-admin-web";
				}
			};

			ManagerConfiguration.init();
			logger.info("create a message endpoint for admin-web");
			endpoint = new JMSMessagingEndpoint(nodeSupport);					

		} catch (MicroarrayException e) {
			if (e.getCause() != null) {
				//The cause has better message, at least when the broker is not available 
				Throwable cause = e.getCause();
				Notification notification = new Notification(cause.getClass().getSimpleName(), cause.getMessage(), Type.ERROR_MESSAGE); 
				notification.show(this.getPage());
			}
			logger.error(e);
			e.printStackTrace();
		} catch (IOException e) {
			logger.error(e);
			e.printStackTrace();
		} catch (IllegalConfigurationException e) {
			logger.error(e);
			e.printStackTrace();
		}
		
		buildMainLayout();
		
		addDetachListener(this);
	}

	@Override
	public void detach(DetachEvent event) {
		// this is executed by default after 15 minutes after the client side
		// has stopped sending heartbeats
		logger.info("clean inactive admin-web session");
		if (endpoint != null) {
			try {
				// delete topics and close the connection
				endpoint.close();
			} catch (JMSException e) {
				logger.error("closing of the messaging endpoint failed", e);
			}
		}
	}

	public MessagingEndpoint getEndpoint() {
		return this.endpoint;
	}
}
