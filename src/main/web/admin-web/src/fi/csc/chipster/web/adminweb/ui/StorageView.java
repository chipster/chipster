package fi.csc.chipster.web.adminweb.ui;

import java.io.IOException;
import java.util.concurrent.locks.Lock;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import com.vaadin.data.Property;
import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.server.AbstractClientConnector;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.Notification;
import com.vaadin.ui.Notification.Type;
import com.vaadin.ui.Panel;
import com.vaadin.ui.ProgressIndicator;

import fi.csc.chipster.web.adminweb.ChipsterAdminUI;
import fi.csc.chipster.web.adminweb.data.StorageAdminAPI;
import fi.csc.chipster.web.adminweb.data.StorageAggregate;
import fi.csc.chipster.web.adminweb.data.StorageAggregateContainer;
import fi.csc.chipster.web.adminweb.data.StorageEntryContainer;
import fi.csc.chipster.web.adminweb.util.StringUtils;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;

@SuppressWarnings("serial")
public class StorageView extends AsynchronousView implements ClickListener, ValueChangeListener, AfterUpdateCallBack {
	
	private static final Logger logger = Logger.getLogger(StorageView.class);

	private static final String DISK_USAGE_BAR_CAPTION = "Total disk usage";
	private StorageEntryTable entryTable;
	private StorageAggregateTable aggregateTable;

	private HorizontalLayout toolbarLayout;

	private StorageEntryContainer entryDataSource;
	private StorageAggregateContainer aggregateDataSource;

	private ChipsterAdminUI app;
	private ProgressIndicator diskUsageBar;
	private HorizontalLayout storagePanels;
	
	private StorageAdminAPI adminEndpoint;

	public StorageView(ChipsterAdminUI app) {

		this.app = app;

		this.addComponent(getToolbar());
		
		this.addComponent(super.getProggressIndicator());
		
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
			
			adminEndpoint = new StorageAdminAPI();
			entryDataSource = new StorageEntryContainer(adminEndpoint);
			aggregateDataSource = new StorageAggregateContainer(adminEndpoint);
			
			entryTable.setContainerDataSource(entryDataSource);
			aggregateTable.setContainerDataSource(aggregateDataSource);
		} catch (JMSException | IOException | IllegalConfigurationException | MicroarrayException | InstantiationException | IllegalAccessException e) {
			logger.error(e);
		}		
	}

	public void update() {
		
		updateStorageAggregates();
		updateStorageTotals();
	}

	/**
	 * Set disk usage. Calls from other threads are allowed.
	 * 
	 * @param usedSpace
	 * @param freeSpace
	 */
	public void setDiskUsage(long usedSpace, long freeSpace) {
		
		//maybe null if the UI thread hasn't initialized this yet
		if (diskUsageBar.getUI() != null) {
			Lock barLock = diskUsageBar.getUI().getSession().getLockInstance();
			barLock.lock();
			try {
				long used = usedSpace;
				long total = usedSpace + freeSpace;
				float division = used / (float)total;
				
				diskUsageBar.setValue(division);
				diskUsageBar.setCaption(DISK_USAGE_BAR_CAPTION + " ( " + 
						StringUtils.getHumanReadable(used) + " / " + StringUtils.getHumanReadable(total) + " )");
				
				if (division > 0.7) {
					diskUsageBar.removeStyleName("ok");
					diskUsageBar.addStyleName("fail");
				} else {
					diskUsageBar.removeStyleName("fail");
					diskUsageBar.addStyleName("ok");
				}
				
				diskUsageBar.markAsDirty();
			} finally {
				barLock.unlock();
			}
		}
	}

	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			toolbarLayout = new HorizontalLayout();
			
			toolbarLayout.addComponent(super.createRefreshButton(this));
			
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

		if (super.isRefreshButton(source)) {
						
			update();
			updateStorageEntries(null);
		}
	}

	public void valueChange(ValueChangeEvent event) {
		Property<?> property = event.getProperty();
		if (property == aggregateTable) {
						
			Object tableValue = aggregateTable.getValue();
			if (tableValue instanceof StorageAggregate) {
				StorageAggregate storageUser = (StorageAggregate) tableValue;
				updateStorageEntries(storageUser.getUsername());
			}			
		}
	}
	
	private void updateStorageEntries(final String username) {
		super.submitUpdate(new Runnable() {

			@Override
			public void run() {				
				entryDataSource.update(StorageView.this, username);
			}			
		}, this);
	}
	
	private void updateStorageAggregates() {
		super.submitUpdate(new Runnable() {

			@Override
			public void run() {				
				aggregateDataSource.update(StorageView.this);
			}			
		}, this);
	}
	
	private void updateStorageTotals() {
		super.submitUpdate(new Runnable() {

			@Override
			public void run() {								
				try {
					Long[] totals = adminEndpoint.getStorageUsage();	
					
					if (totals != null) {
						StorageView.this.setDiskUsage(totals[0], totals[1]);
					} else {
						Notification.show("Timeout", "Chipster filebroker server doesn't respond", Type.ERROR_MESSAGE);
						logger.error("timeout while waiting storage usage totals");
					}

				} catch (JMSException | InterruptedException e) {
					logger.error(e);
				}			
			}			
		}, this);
	}
	
	public ChipsterAdminUI getApp() {
		return app;
	}

	public void delete(Object itemId) {
		//TODO are you sure?
		try {
			adminEndpoint.deleteRemoteSession(entryDataSource.getItem(itemId).getBean().getID());
		} catch (JMSException e) {
			logger.warn("could not delete session", e);
			return;
		}
		
		entryDataSource.removeItem(itemId);
		
		aggregateDataSource.update(this);
		updateStorageTotals();
	}
	
	/**
	 * Calling from background threads allowed
	 */
	@Override
	public void updateDone() {
					
				
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
				aggregateTable.setVisibleColumns(StorageAggregateContainer.NATURAL_COL_ORDER);
				aggregateTable.setColumnHeaders(StorageAggregateContainer.COL_HEADERS_ENGLISH);				

			} finally {
				aggregateTableLock.unlock();
			}
		}
	}

	public AbstractClientConnector getEntryTable() {
		return entryTable;
	}

	public void clean() {
		if (adminEndpoint != null) {
			adminEndpoint.clean();
		}
	}
}
