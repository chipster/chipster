package fi.csc.microarray.manager.web.ui;

import com.vaadin.data.Property;
import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.terminal.ThemeResource;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.HorizontalLayout;

import fi.csc.microarray.manager.web.ChipsterAdminApplication;
import fi.csc.microarray.manager.web.data.StorageAggregate;
import fi.csc.microarray.manager.web.data.StorageAggregateContainer;
import fi.csc.microarray.manager.web.data.StorageEntryContainer;

public class StorageView extends HorizontalLayout implements ClickListener, ValueChangeListener {

	private StorageEntryTable entryTable;
	private StorageAggregateTable aggregateTable;

	private HorizontalLayout toolbarLayout;

	private Button refreshButton = new Button("Refresh");

	private StorageEntryContainer entryDataSource;
	private StorageAggregateContainer aggregateDataSource;

	private ChipsterAdminApplication app;


	public StorageView(ChipsterAdminApplication app) {

		this.app = app;

		this.addComponent(getToolbar());

		entryTable = new StorageEntryTable(this);
		aggregateTable = new StorageAggregateTable(this);

		this.addComponent(aggregateTable);
		this.addComponent(entryTable);
		this.setExpandRatio(entryTable, 1);
		this.setSizeFull();
	}

	public void loadData() throws InstantiationException, IllegalAccessException {

		entryDataSource = new StorageEntryContainer();
		aggregateDataSource = new StorageAggregateContainer(entryDataSource);

		entryTable.setContainerDataSource(entryDataSource);
		aggregateTable.setContainerDataSource(aggregateDataSource);

		entryDataSource.update(this);
		StorageAggregate totalItem = aggregateDataSource.update(this);
		aggregateTable.select(totalItem);
		
		entryTable.setVisibleColumns(StorageEntryContainer.NATURAL_COL_ORDER);
		entryTable.setColumnHeaders(StorageEntryContainer.COL_HEADERS_ENGLISH);

		aggregateTable.setVisibleColumns(StorageAggregateContainer.NATURAL_COL_ORDER);
		aggregateTable.setColumnHeaders(StorageAggregateContainer.COL_HEADERS_ENGLISH);
	}

	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			toolbarLayout = new HorizontalLayout();
			refreshButton.addListener((ClickListener)this);
			toolbarLayout.addComponent(refreshButton);

			refreshButton.setIcon(new ThemeResource("../runo/icons/32/reload.png"));

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
		}
	}

	public void valueChange(ValueChangeEvent event) {
		Property property = event.getProperty();
		if (property == aggregateTable) {
			StorageAggregate selection = (StorageAggregate) aggregateTable.getValue();

			if (aggregateDataSource.TOTAL_USERNAME.equals(selection.getUsername())) {
				entryDataSource.showUser(null);
			} else {
				entryDataSource.showUser(selection.getUsername());
			}
		}
	}

	public ChipsterAdminApplication getApp() {
		return app;
	}
}
