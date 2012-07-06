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

import fi.csc.microarray.manager.web.ChipsterAdminApplication;
import fi.csc.microarray.manager.web.data.JobsContainer;

public class JobsView extends VerticalLayout implements ClickListener, ValueChangeListener  {
	
	private HorizontalLayout toolbarLayout;

	private Button refreshButton = new Button("Refresh");

	private JobsTable table;
	private JobsContainer dataSource;

	public JobsView(ChipsterAdminApplication app) {
		
		dataSource = new JobsContainer(); 
				
		table = new JobsTable(this);
		table.setContainerDataSource(dataSource);

		table.setVisibleColumns(JobsContainer.NATURAL_COL_ORDER);
		table.setColumnHeaders(JobsContainer.COL_HEADERS_ENGLISH);
		
		this.addComponent(getToolbar());
		this.addComponent(table);

		setSizeFull();
		this.setExpandRatio(table, 1);
	}

	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			
			toolbarLayout = new HorizontalLayout();
			
			Label spaceEater = new Label(" ");
			toolbarLayout.addComponent(spaceEater);
			toolbarLayout.setExpandRatio(spaceEater, 1);
			
			refreshButton.addListener((ClickListener)this);
			refreshButton.setIcon(new ThemeResource("../runo/icons/32/reload.png"));
			toolbarLayout.addComponent(refreshButton);
			
			toolbarLayout.setWidth("100%");
			toolbarLayout.setStyleName("toolbar");
		}

		return toolbarLayout;
	}

	public void buttonClick(ClickEvent event) {
		final Button source = event.getButton();

		if (source == refreshButton) {
			dataSource.update();
		} 
	}

	public void valueChange(ValueChangeEvent event) {
		Property property = event.getProperty();
		if (property == table) {
			//			Item item = personList.getItem(personList.getValue());
			//			if (item != personForm.getItemDataSource()) {
			//				personForm.setItemDataSource(item);
			//			}
		}
	}

	public void cancel(Object itemId) {
		dataSource.removeItem(itemId);
	}
}
