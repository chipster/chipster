package fi.csc.microarray.manager.web.ui;

import com.vaadin.data.Property;
import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.terminal.ThemeResource;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.Panel;
import com.vaadin.ui.ProgressIndicator;
import com.vaadin.ui.VerticalLayout;

import fi.csc.microarray.manager.web.ChipsterAdminApplication;
import fi.csc.microarray.manager.web.data.StorageAggregate;
import fi.csc.microarray.manager.web.data.StorageAggregateContainer;
import fi.csc.microarray.manager.web.data.StorageEntryContainer;
import fi.csc.microarray.manager.web.util.StringUtils;

public class StorageView extends VerticalLayout implements ClickListener, ValueChangeListener {

	private static final String DISK_USAGE_BAR_CAPTION = "Total disk usage";
	private StorageEntryTable entryTable;
	private StorageAggregateTable aggregateTable;

	private HorizontalLayout toolbarLayout;

	private Button refreshButton = new Button("Refresh");

	private StorageEntryContainer entryDataSource;
	private StorageAggregateContainer aggregateDataSource;

	private ChipsterAdminApplication app;
	private ProgressIndicator diskUsageBar;
	private HorizontalLayout storagePanels;


	public StorageView(ChipsterAdminApplication app) {

		this.app = app;

		this.addComponent(getToolbar());

		
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
		
		aggregatePanel.setWidth(300, UNITS_PIXELS);
		aggregatePanel.setHeight(100, UNITS_PERCENTAGE);
		entryPanel.setSizeFull();
		
		aggregatePanel.setContent(aggregatePanelLayout);
		entryPanel.setContent(entryPanelLayout);
		
		storagePanels = new HorizontalLayout();
		storagePanels.setSizeFull();
					
		storagePanels.addComponent(aggregatePanel);
		storagePanels.addComponent(entryPanel);
		
		//storagePanels.setExpandRatio(aggregatePanel, 1);
		storagePanels.setExpandRatio(entryPanel, 1);
		
		this.setSizeFull();
		this.addComponent(storagePanels);		
		this.setExpandRatio(storagePanels, 1);
	}

	public void loadData() throws InstantiationException, IllegalAccessException {

		entryDataSource = new StorageEntryContainer();
		aggregateDataSource = new StorageAggregateContainer(entryDataSource);

		entryTable.setContainerDataSource(entryDataSource);
		aggregateTable.setContainerDataSource(aggregateDataSource);

		entryDataSource.update(this);
		
		StorageAggregate totalItem = aggregateDataSource.update(this);
		aggregateTable.select(totalItem);
		
		updateDiskUsageBar();
		
		entryTable.setVisibleColumns(StorageEntryContainer.NATURAL_COL_ORDER);
		entryTable.setColumnHeaders(StorageEntryContainer.COL_HEADERS_ENGLISH);

		aggregateTable.setVisibleColumns(StorageAggregateContainer.NATURAL_COL_ORDER);
		aggregateTable.setColumnHeaders(StorageAggregateContainer.COL_HEADERS_ENGLISH);
	}

	private void updateDiskUsageBar() {
		
		long used = aggregateDataSource.getDiskUsage();
		long total = aggregateDataSource.getDiskAvailable() + used;
		
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
	}

	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			toolbarLayout = new HorizontalLayout();
			refreshButton.addListener((ClickListener)this);
			toolbarLayout.addComponent(refreshButton);

			refreshButton.setIcon(new ThemeResource("../runo/icons/32/reload.png"));
			
			Label spaceEater = new Label(" ");
			toolbarLayout.addComponent(spaceEater);
			toolbarLayout.setExpandRatio(spaceEater, 1);
			
			diskUsageBar = new ProgressIndicator(0f);
			diskUsageBar.setCaption(DISK_USAGE_BAR_CAPTION);
			diskUsageBar.setStyleName("big");
			
			diskUsageBar.setWidth(300, UNITS_PIXELS);
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
			entryDataSource.update(this);
			aggregateDataSource.update(this);
			updateDiskUsageBar();
		}
	}

	public void valueChange(ValueChangeEvent event) {
		Property property = event.getProperty();
		if (property == aggregateTable) {
			updateEntryFilter();
			StorageAggregate selection = (StorageAggregate) aggregateTable.getValue();

			if (aggregateDataSource.TOTAL_USERNAME.equals(selection.getUsername())) {
				entryDataSource.showUser(null);
			} else {
			}
		}
	}
	
	public void updateEntryFilter() {
		
		StorageAggregate selection = (StorageAggregate) aggregateTable.getValue();

		if (aggregateDataSource.TOTAL_USERNAME.equals(selection.getUsername())) {
			entryDataSource.showUser(null);
		} else {
			entryDataSource.showUser(selection.getUsername());
		}
	}
	
	public ChipsterAdminApplication getApp() {
		return app;
	}

	public void delete(Object itemId) {
		//TODO are you sure?
		//TODO remove from the server
		entryDataSource.removeItem(itemId);
		
		aggregateDataSource.update(this);
		updateEntryFilter();
		updateDiskUsageBar();
	}
}
