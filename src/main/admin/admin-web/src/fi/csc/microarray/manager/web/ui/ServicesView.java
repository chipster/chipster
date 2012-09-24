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
import com.vaadin.ui.VerticalLayout;
import com.vaadin.ui.Window.Notification;

import fi.csc.microarray.manager.web.ChipsterAdminApplication;
import fi.csc.microarray.manager.web.data.ServiceContainer;

public class ServicesView extends VerticalLayout implements ClickListener, ValueChangeListener {
	
	private ServicesTable table;
	private HorizontalLayout toolbarLayout;
	
 	private Button refreshButton = new Button("Refresh");

	private ServiceContainer dataSource;
	private ChipsterAdminApplication app;


	public ServicesView(ChipsterAdminApplication app) {
		
		this.app = app;
		
		this.addComponent(getToolbar());

		table = new ServicesTable(this);

		this.addComponent(table);
		this.setExpandRatio(table, 1);
		this.setSizeFull();
	}
	
	public void loadData() throws InstantiationException, IllegalAccessException {
		dataSource = new ServiceContainer();
		table.setContainerDataSource(dataSource);
		dataSource.update(this);
	}
	
	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			toolbarLayout = new HorizontalLayout();
			
			refreshButton.setIcon(new ThemeResource("../runo/icons/32/reload.png"));
			refreshButton.addListener((ClickListener)this);
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
			dataSource.update(this);
		}
	}

	public void valueChange(ValueChangeEvent event) {
		Property property = event.getProperty();
		if (property == table) {
			//Nothing to do yet
		}
	}

	public ChipsterAdminApplication getApp() {
		return app;
	}

	public void dataUpdated() {
		table.setVisibleColumns(ServiceContainer.NATURAL_COL_ORDER);
		table.setColumnHeaders(ServiceContainer.COL_HEADERS_ENGLISH);
		
		getApp().getMainWindow().showNotification(
				"Found "
						+ table.getContainerDataSource().size() + " nodes",
						Notification.TYPE_TRAY_NOTIFICATION);
	}
}
