package fi.csc.chipster.web.adminweb.ui;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.locks.Lock;

import com.vaadin.data.Property;
import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.server.AbstractClientConnector;
import com.vaadin.server.ThemeResource;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.Panel;
import com.vaadin.ui.ProgressIndicator;
import com.vaadin.ui.VerticalLayout;

import fi.csc.chipster.web.adminweb.ChipsterAdminUI;
import fi.csc.chipster.web.adminweb.data.StorageAggregate;
import fi.csc.chipster.web.adminweb.data.StorageAggregateContainer;
import fi.csc.chipster.web.adminweb.data.StorageEntryContainer;
import fi.csc.chipster.web.adminweb.util.StringUtils;

public class StorageView extends VerticalLayout implements ClickListener, ValueChangeListener {

	private static final String DISK_USAGE_BAR_CAPTION = "Total disk usage";
	private StorageEntryTable entryTable;
	private StorageAggregateTable aggregateTable;

	private HorizontalLayout toolbarLayout;

	private Button refreshButton = new Button("Refresh");

	private StorageEntryContainer entryDataSource;
	private StorageAggregateContainer aggregateDataSource;

	private ChipsterAdminUI app;
	private ProgressIndicator diskUsageBar;
	private HorizontalLayout storagePanels;
	
	private ProgressIndicator progressIndicator = new ProgressIndicator(0.0f);
	private boolean updateDone;
	private static final int POLLING_INTERVAL = 100;

	public StorageView(ChipsterAdminUI app) {

		this.app = app;

		this.addComponent(getToolbar());
		
		progressIndicator.setWidth(100, Unit.PERCENTAGE);
		this.addComponent(progressIndicator);
		
		entryTable = new StorageEntryTable(this);
		aggregateTable = new StorageAggregateTable(this);

		
		HorizontalLayout aggregatePanelLayout = new HorizontalLayout();
		HorizontalLayout entryPanelLayout = new HorizontalLayout();
		
		aggregatePanelLayout.setSizeFull();
		entryPanelLayout.setSizeFull();
		
		aggregatePanelLayout.addComponent(aggregateTable);
		entryPanelLayout.addComponent(entryTable);
		
		aggregatePanelLayout.setExpandRatio(aggregateTable, 1);
		entryPanelLayout.setExpandRatio(entryTable, 1);
		
		Panel aggregatePanel = new Panel("Disk usage by user");
		Panel entryPanel = new Panel("Stored sessions");
		
		aggregatePanel.setWidth(300, Unit.PIXELS);
		aggregatePanel.setHeight(100, Unit.PERCENTAGE);
		entryPanel.setSizeFull();
		
		aggregatePanel.setContent(aggregatePanelLayout);
		entryPanel.setContent(entryPanelLayout);
		
		storagePanels = new HorizontalLayout();
		storagePanels.setSizeFull();
					
		storagePanels.addComponent(aggregatePanel);
		storagePanels.addComponent(entryPanel);
		
		storagePanels.setExpandRatio(entryPanel, 1);
		
		this.setSizeFull();
		this.addComponent(storagePanels);		
		this.setExpandRatio(storagePanels, 1);
		
		try {
			entryDataSource = new StorageEntryContainer();
			aggregateDataSource = new StorageAggregateContainer(entryDataSource);
			
			entryTable.setContainerDataSource(entryDataSource);
			aggregateTable.setContainerDataSource(aggregateDataSource);
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}		
	}

	public void update() {
		
		aggregateDataSource.update(this);
		waitForUpdate();
	}
	
	private void waitForUpdate() {
		
		//Disable during data update avoid concurrent modification
		refreshButton.setEnabled(false);
		
		updateDone = false;
		
		progressIndicator.setPollingInterval(POLLING_INTERVAL);
		
		ExecutorService execService = Executors.newCachedThreadPool();
		execService.execute(new Runnable() {
			public void run() {
				
				try {
					/* Separate delay from what happens in the Container, because communication between
					 * threads is messy. Nevertheless, these delays should have approximately same duration
					 * to prevent user from starting several background updates causing concurrent modifications.   
					 */
					final int DELAY = 300; 				
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

	private void updateDiskUsageBar() {
		
//		long used = aggregateDataSource.getDiskUsage();
//		long total = aggregateDataSource.getDiskAvailable() + used;
		
		long used = 250000000000l;
		long total = 500000000000l;
		
		diskUsageBar.setValue(used / (float)total);
		diskUsageBar.setCaption(DISK_USAGE_BAR_CAPTION + " ( " + 
				StringUtils.getHumanReadable(used) + " / " + StringUtils.getHumanReadable(total) + " )");
		
		if (used / (float)total > 0.7) {
			diskUsageBar.removeStyleName("ok");
			diskUsageBar.addStyleName("fail");
		} else {
			diskUsageBar.removeStyleName("fail");
			diskUsageBar.addStyleName("ok");
		}
		
		diskUsageBar.markAsDirty();
	}

	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			toolbarLayout = new HorizontalLayout();
			refreshButton.addClickListener((ClickListener)this);
			toolbarLayout.addComponent(refreshButton);

			refreshButton.setIcon(new ThemeResource("../runo/icons/32/reload.png"));
			
			Label spaceEater = new Label(" ");
			toolbarLayout.addComponent(spaceEater);
			toolbarLayout.setExpandRatio(spaceEater, 1);
			
			diskUsageBar = new ProgressIndicator(0f);
			diskUsageBar.setCaption(DISK_USAGE_BAR_CAPTION);
			diskUsageBar.setStyleName("big");
			diskUsageBar.setPollingInterval(Integer.MAX_VALUE);
			
			diskUsageBar.setWidth(300, Unit.PIXELS);
			toolbarLayout.addComponent(diskUsageBar);
			toolbarLayout.setExpandRatio(diskUsageBar, 1);
					
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
//			entryDataSource.update(this);
//			aggregateDataSource.update(this);
//			updateDiskUsageBar();
		}
	}

	public void valueChange(ValueChangeEvent event) {
		Property<?> property = event.getProperty();
		if (property == aggregateTable) {
						
			Object tableValue = aggregateTable.getValue();
			if (tableValue instanceof StorageAggregate) {
				StorageAggregate storageUser = (StorageAggregate) tableValue;
				
				entryDataSource.update(this, storageUser.getUsername());
				waitForUpdate();
			}			
		}
	}
	
	public ChipsterAdminUI getApp() {
		return app;
	}

	public void delete(Object itemId) {
		//TODO are you sure?
		//TODO remove from the server
		entryDataSource.removeItem(itemId);
		
		aggregateDataSource.update(this);
		updateDiskUsageBar();
	}
	
	/**
	 * Calling from background threads allowed
	 */
	public void entryUpdateDone() {
					
				
		if (entryTable.getUI() != null) {
			Lock entryTableLock = entryTable.getUI().getSession().getLockInstance();
			entryTableLock.lock();
			try {

				entryTable.setVisibleColumns(StorageEntryContainer.NATURAL_COL_ORDER);
				entryTable.setColumnHeaders(StorageEntryContainer.COL_HEADERS_ENGLISH);

			} finally {
				entryTableLock.unlock();
			}
		}
		
		if (aggregateTable.getUI() != null) {
			Lock aggregateTableLock = aggregateTable.getUI().getSession().getLockInstance();
			aggregateTableLock.lock();
			try {						
				//StorageAggregate totalItem = aggregateDataSource.update(this);

				//aggregateTable.select(totalItem);
				aggregateTable.setVisibleColumns(StorageAggregateContainer.NATURAL_COL_ORDER);
				aggregateTable.setColumnHeaders(StorageAggregateContainer.COL_HEADERS_ENGLISH);				

			} finally {
				aggregateTableLock.unlock();
			}
		}
		
		if (progressIndicator.getUI() != null) {
			Lock proggressIndicatorLock = progressIndicator.getUI().getSession().getLockInstance();
			proggressIndicatorLock.lock();
			try {						
				updateDiskUsageBar();
			} finally {
				proggressIndicatorLock.unlock();
			}
		}

		this.updateDone = true;
	}

	public AbstractClientConnector getEntryTable() {
		return entryTable;
	}
}
