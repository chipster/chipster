package fi.csc.microarray.manager.web.ui;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

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

import fi.csc.microarray.manager.web.ChipsterAdminUI;
import fi.csc.microarray.manager.web.data.ServiceContainer;

public class ServicesView extends VerticalLayout implements ClickListener, ValueChangeListener {

	private ServicesTable table;
	private HorizontalLayout toolbarLayout;

	private Button refreshButton = new Button("Refresh");

	private ServiceContainer dataSource;
	private ChipsterAdminUI app;

	private ProgressIndicator progressIndicator = new ProgressIndicator(0.0f);
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
		
		dataSource.update(this);
		
		progressIndicator.setPollingInterval(POLLING_INTERVAL);
		
		ExecutorService execService = Executors.newCachedThreadPool();
		execService.execute(new Runnable() {
			public void run() {
				
				try {
					/* Separate delay from what happens in the ServiceContainer, because communication between
					 * threads is messy. Nevertheless, these delays should have approximately same duration
					 * to prevent user from starting several background updates causing concurrent modifications.   
					 */
					final int DELAY = 50; 				
					for (int i = 0; i <= DELAY; i++) {

						//First case happens in initialisation, second if another view is chosen during the data update 
						if (progressIndicator.getUI() != null && progressIndicator.getUI().getSession().getLock() != null ) {
							
							//Component has to be locked before modification from background thread
							progressIndicator.getUI().getSession().getLock().lock();					
							try {
								progressIndicator.setValue((float)i/DELAY);
							} finally {
								progressIndicator.getUI().getSession().getLock().unlock();
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
					progressIndicator.setPollingInterval(Integer.MAX_VALUE);
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
}
