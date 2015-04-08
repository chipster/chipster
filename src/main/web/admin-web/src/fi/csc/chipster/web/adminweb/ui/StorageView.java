package fi.csc.chipster.web.adminweb.ui;

import java.io.IOException;

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
import com.vaadin.ui.ProgressBar;

import fi.csc.chipster.web.adminweb.ChipsterAdminUI;
import fi.csc.chipster.web.adminweb.data.StorageAggregateContainer;
import fi.csc.chipster.web.adminweb.data.StorageEntryContainer;
import fi.csc.chipster.web.adminweb.util.Notificationutil;
import fi.csc.chipster.web.adminweb.util.StringUtils;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.admin.StorageAdminAPI;
import fi.csc.microarray.messaging.admin.StorageAggregate;

@SuppressWarnings("serial")
public class StorageView extends AsynchronousView implements ClickListener, ValueChangeListener, AfterUpdateCallBack {
	
	private static final Logger logger = Logger.getLogger(StorageView.class);

	private static final String DISK_USAGE_BAR_CAPTION = "Total disk usage";
	private StorageEntryTable entryTable;
	private StorageAggregateTable aggregateTable;

	private HorizontalLayout toolbarLayout;

	private StorageEntryContainer entryDataSource;
	private StorageAggregateContainer aggregateDataSource;

	private ProgressBar diskUsageBar;
	private HorizontalLayout storagePanels;
	
	private StorageAdminAPI adminEndpoint;

	public StorageView(ChipsterAdminUI app) {
		super(app);

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
			
			adminEndpoint = new StorageAdminAPI(app.getEndpoint());
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
	public void setDiskUsage(final long usedSpace, final long freeSpace) {

		this.updateUI(new Runnable() {
			public void run() {
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
			}
		});
	}

	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			toolbarLayout = new HorizontalLayout();
			
			toolbarLayout.addComponent(super.createRefreshButton(this));
			
			Label spaceEater = new Label(" ");
			toolbarLayout.addComponent(spaceEater);
			toolbarLayout.setExpandRatio(spaceEater, 1);
			
			diskUsageBar = new ProgressBar(0f);
			diskUsageBar.setCaption(DISK_USAGE_BAR_CAPTION);
			diskUsageBar.setStyleName("big");
			
			diskUsageBar.setWidth(300, Unit.PIXELS);
			toolbarLayout.addComponent(diskUsageBar);
			toolbarLayout.setExpandRatio(diskUsageBar, 1);
					
			toolbarLayout.addComponent(super.getApp().getTitle());	

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

	public void delete(Object itemId) {
		//TODO are you sure?
		try {
			adminEndpoint.deleteRemoteSession(entryDataSource.getItem(itemId).getBean().getID());
		} catch (JMSException | MicroarrayException e) {			
			Notificationutil.showFailNotification(e.getClass().getSimpleName(), e.getMessage());
			logger.warn("could not delete session", e);
			return;
		}
		
		entryDataSource.removeItem(itemId);
		
		aggregateDataSource.update(this);
		updateStorageTotals();
	}
	
	@Override
	public void updateDone() {
				
		entryTable.setVisibleColumns(StorageEntryContainer.NATURAL_COL_ORDER);
		entryTable.setColumnHeaders(StorageEntryContainer.COL_HEADERS_ENGLISH);

		aggregateTable.setVisibleColumns(StorageAggregateContainer.NATURAL_COL_ORDER);
		aggregateTable.setColumnHeaders(StorageAggregateContainer.COL_HEADERS_ENGLISH);	
	}

	public AbstractClientConnector getEntryTable() {
		return entryTable;
	}
}
