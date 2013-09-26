package fi.csc.chipster.web.adminweb.ui;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.locks.Lock;

import com.vaadin.data.Property;
import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.server.ThemeResource;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.ProgressIndicator;
import com.vaadin.ui.VerticalLayout;

import fi.csc.chipster.web.adminweb.ChipsterAdminUI;
import fi.csc.chipster.web.adminweb.data.ServiceContainer;
import fi.csc.chipster.web.adminweb.data.ServiceEntry;
import fi.csc.microarray.messaging.AdminAPI.NodeStatus.Status;

public class ServicesView extends VerticalLayout implements ClickListener, ValueChangeListener {

	private ServicesTable table;
	private HorizontalLayout toolbarLayout;

	private Button refreshButton = new Button("Refresh");

	private ServiceContainer dataSource;
	private ChipsterAdminUI app;

	private ProgressIndicator progressIndicator = new ProgressIndicator(0.0f);
	private boolean updateDone;
	private static final int POLLING_INTERVAL = 100; 


	public ServicesView(ChipsterAdminUI app) {

		this.app = app;

		this.addComponent(getToolbar());

		progressIndicator.setWidth(100, Unit.PERCENTAGE);
		this.addComponent(progressIndicator);

		table = new ServicesTable(this);

		this.addComponent(table);
		this.setExpandRatio(table, 1);
		this.setSizeFull();

		try {
			dataSource = new ServiceContainer();
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}		
		table.setContainerDataSource(dataSource);			

		table.setVisibleColumns(ServiceContainer.NATURAL_COL_ORDER);
		table.setColumnHeaders(ServiceContainer.COL_HEADERS_ENGLISH);
		
		table.setSortContainerPropertyId(ServiceContainer.NAME);
	}

	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			toolbarLayout = new HorizontalLayout();

			refreshButton.setIcon(new ThemeResource("../runo/icons/32/reload.png"));
			refreshButton.addClickListener((ClickListener)this);
			toolbarLayout.addComponent(refreshButton);

			Label spaceEater = new Label(" ");
			toolbarLayout.addComponent(spaceEater);
			toolbarLayout.setExpandRatio(spaceEater, 1);

			toolbarLayout.addComponent(app.getTitle());	

			toolbarLayout.setWidth("100%");
			toolbarLayout.setStyleName("toolbar");
		}
		return toolbarLayout;

	}

	public void buttonClick(ClickEvent event) {
		final Button source = event.getButton();

		if (source == refreshButton) {
			update();
		}
	}

	public void update() {
		//Disable during data update avoid concurrent modification
		refreshButton.setEnabled(false);
		
		updateDone = false;
		dataSource.update(this);
		
		progressIndicator.setPollingInterval(POLLING_INTERVAL);
		
		ExecutorService execService = Executors.newCachedThreadPool();
		execService.execute(new Runnable() {
			public void run() {
				
				try {
					/* Separate delay from what happens in the ServiceContainer, because communication between
					 * threads is messy. Nevertheless, these delays should have approximately same duration
					 * to prevent user from starting several background updates causing concurrent modifications.
					 * 
					 * When the broker is not available, the connection will timeout in 30 seconds. When this connection
					 * error happens, a notification is created, but there isn't simple way for the server to push it to
					 * client side. ProgressIndicator however, creates page updates and shows the notification if this 
					 * timeout is longer than the connection timeout.
					 */
					final int DELAY = 350; 				
					for (int i = 0; i <= DELAY; i++) {
						
						if (updateDone) {							
							break;
						}

						//This happens in initialisation 
						if (progressIndicator.getUI() != null ) {
							
							Lock indicatorLock = progressIndicator.getUI().getSession().getLockInstance();
							
							//Component has to be locked before modification from background thread
							indicatorLock.lock();					
							try {
								progressIndicator.setValue((float)i/DELAY);
							} finally {
								indicatorLock.unlock();
							}
						}

						try {
							Thread.sleep(100);
						} catch (InterruptedException e) {
							//Just continue
						}
					}
					
				} finally {
					refreshButton.setEnabled(true);
					
					if (progressIndicator.getUI() != null) {
						Lock indicatorLock = progressIndicator.getUI().getSession().getLockInstance();

						indicatorLock.lock();					
						try {
							progressIndicator.setValue(1.0f);
							progressIndicator.setPollingInterval(Integer.MAX_VALUE);
						} finally {
							indicatorLock.unlock();
						}
					}
				}
			}
		});
	}

	public void valueChange(ValueChangeEvent event) {
		Property<?> property = event.getProperty();
		if (property == table) {
			//Nothing to do yet
		}
	}

	public ChipsterAdminUI getApp() {
		return app;
	}

	public ServicesTable getTable() {
		return table;
	}

	
	/**
	 * Calling from background threads allowed
	 */
	public void updateDone() {
			
		if (table.getUI() != null) {
			Lock tableLock = table.getUI().getSession().getLockInstance();
			tableLock.lock();
			try {

				for (ServiceEntry entry : dataSource.getItemIds()) {
					if (Status.UNKNOWN.equals(entry.getStatus())) {
						entry.setStatus(Status.DOWN);
					}
				}

				table.markAsDirtyRecursive();						

			} finally {
				tableLock.unlock();
			}
		}

		this.updateDone = true;
	}
}
