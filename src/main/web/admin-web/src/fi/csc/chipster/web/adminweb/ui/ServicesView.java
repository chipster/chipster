package fi.csc.chipster.web.adminweb.ui;

import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;

import fi.csc.chipster.web.adminweb.ChipsterAdminUI;
import fi.csc.chipster.web.adminweb.data.ServiceContainer;
import fi.csc.chipster.web.adminweb.data.ServiceEntry;
import fi.csc.microarray.messaging.admin.AdminAPI.NodeStatus.Status;

public class ServicesView extends AsynchronousView implements ClickListener, AfterUpdateCallBack {

	private ServicesTable table;
	private HorizontalLayout toolbarLayout;

	private ServiceContainer dataSource;

	public ServicesView(ChipsterAdminUI app) {
		
		super(app, 5);

		this.addComponent(getToolbar());

		this.addComponent(getProggressIndicator());

		table = new ServicesTable(this);

		this.addComponent(table);
		this.setExpandRatio(table, 1);
		this.setSizeFull();

		try {
			dataSource = new ServiceContainer(this);
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
			
			toolbarLayout.addComponent(super.createRefreshButton(this));

			Label spaceEater = new Label(" ");
			toolbarLayout.addComponent(spaceEater);
			toolbarLayout.setExpandRatio(spaceEater, 1);

			toolbarLayout.addComponent(getApp().getTitle());	

			toolbarLayout.setWidth("100%");
			toolbarLayout.setStyleName("toolbar");
		}
		return toolbarLayout;

	}

	public void buttonClick(ClickEvent event) {
		final Button source = event.getButton();

		if (isRefreshButton(source)) {
			update();
		}
	}

	public void update() {
		
		super.submitUpdate(new Runnable() {

			@Override
			public void run() {				
				dataSource.update();
			}			
		}, this);				
	}

	public ServicesTable getTable() {
		return table;
	}
	
	@Override
	public void updateDone() {
		for (ServiceEntry entry : dataSource.getItemIds()) {
			if (Status.UNKNOWN.equals(entry.getStatus())) {
				entry.setStatus(Status.DOWN);
			}
		}
		table.markAsDirtyRecursive();						
	}
}
