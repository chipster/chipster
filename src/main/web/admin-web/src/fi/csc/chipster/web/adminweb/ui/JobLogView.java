package fi.csc.chipster.web.adminweb.ui;

import java.util.LinkedList;

import org.apache.log4j.Logger;
import org.hibernate.exception.GenericJDBCException;

import com.vaadin.data.Property;
import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.server.ThemeResource;
import com.vaadin.shared.ui.label.ContentMode;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.CheckBox;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.Notification;
import com.vaadin.ui.VerticalLayout;
import com.vaadin.ui.Window;

import fi.csc.chipster.web.adminweb.ChipsterAdminUI;
import fi.csc.chipster.web.adminweb.data.DateContainerFilter;
import fi.csc.chipster.web.adminweb.data.JobLogContainer;

public class JobLogView extends VerticalLayout implements ClickListener  {
	
	private static final Logger logger = Logger.getLogger(JobLogView.class);
	
	private HorizontalLayout toolbarLayout;

	private Button addFilterButton = new Button();

	private JobLogTable table;
	private JobLogContainer dataSource;

	private ChipsterAdminUI app;

	private LinkedList<JobLogFilter> filters;


	private HorizontalLayout filterLayout;

	private CheckBox ignoreTestAccounts;


	public JobLogView(ChipsterAdminUI app) {

		this.app = app;
		// do this before data source is attached to avoid one data update
		setSizeFull();
		table = new JobLogTable(this);
		
		this.addComponent(getToolbar());
		this.addComponent(table);

		this.setExpandRatio(table, 1);
		
		try {
			dataSource = new JobLogContainer(this);			
			table.setContainerDataSource(dataSource);					

			table.setVisibleColumns(JobLogContainer.NATURAL_COL_ORDER);
			table.setColumnHeaders(JobLogContainer.COL_HEADERS_ENGLISH);

			table.setSortAscending(false);
			table.setSortContainerPropertyId(JobLogContainer.END_TIME);
			
			addFilter(JobLogContainer.END_TIME, DateContainerFilter.getToday());
			
		} catch (GenericJDBCException e) {
			logger.error("unable to read job database", e);
			//FIXME Show exception message and hide or disable all database based content 
			return;
		}		
	}
	
	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			
			toolbarLayout = new HorizontalLayout();
			filterLayout = new HorizontalLayout();
			
			toolbarLayout.addComponent(filterLayout);
			addFilterButton.addClickListener((ClickListener)this);
			addFilterButton.setIcon(new ThemeResource("crystal/edit_add.png"));
			addFilterButton.setDescription("Add another filter");
			addFilterButton.addStyleName("search-button");
			
			HorizontalLayout buttonBorder = new HorizontalLayout();
			buttonBorder.addStyleName("search-filter-bg");
			buttonBorder.addComponent(addFilterButton);
			toolbarLayout.addComponent(buttonBorder);
			
			Button searchButton = new Button();
			searchButton.setIcon(new ThemeResource("crystal/mail_find.png"));
			searchButton.setDescription("Search");
			toolbarLayout.addComponent(searchButton);
			
			searchButton.addClickListener(new Button.ClickListener() {

				public void buttonClick(ClickEvent event) {

					update();
				}
			});
			
			/*
			 * Current way of filtering
			 * strings is slow, because H2 doesn't use index for these SQL
			 * queries (WHERE NOT username = '').
			 */
			ignoreTestAccounts = new CheckBox("Ignore test accounts", false);
			ignoreTestAccounts.addStyleName("toolbar-component");
			toolbarLayout.addComponent(ignoreTestAccounts);
						
			ignoreTestAccounts.addValueChangeListener(new ValueChangeListener() {

				@Override
				public void valueChange(ValueChangeEvent arg0) {
					update();
				}
			});
			
			Label spaceEater = new Label(" ");
			toolbarLayout.addComponent(spaceEater);
			toolbarLayout.setExpandRatio(spaceEater, 1);

			toolbarLayout.addComponent(app.getTitle());	
			
			toolbarLayout.setWidth("100%");
			toolbarLayout.setStyleName("toolbar");
		}

		return toolbarLayout;
	}

	private void addFilter(String column, String value) {
		if (filters == null) {
			filters = new LinkedList<JobLogFilter>();
		}
		JobLogFilter filter = new JobLogFilter(this, column, value);
		filters.add(filter);
		filterLayout.addComponent(filter);
	}

	public void buttonClick(ClickEvent event) {
		final Button source = event.getButton();

		if (source == addFilterButton) {
			addFilter(null, null);
		}
	}

	public void update() {

		updateContainerFilters();

		Notification.show("Found "
						+ table.getContainerDataSource().size() + " item(s)", Notification.Type.TRAY_NOTIFICATION);
	}
	
	private void updateContainerFilters() {
		dataSource.removeAllContainerFilters();
		dataSource.setIgnoreTestAccounts(ignoreTestAccounts.getValue());

		for (JobLogFilter iteratedSearch : filters) {
			if (iteratedSearch.getContainerFilter() != null) {
				dataSource.addContainerFilter(iteratedSearch.getContainerFilter());
			}
		}
		
		table.refreshRowCache();
	}

	public void clearFilter(JobLogFilter filter) {

		if (filters.size() > 1) {
			filterLayout.removeComponent(filter);
			filters.remove(filter);
		} else {
			filters.get(0).clear();
		}

		//It's not possible to remove just one filter, so let's do it in the hard way
		updateContainerFilters();
	}

	public void showOutput(Object itemId) {
		
		String output = "";
		Property<?> outputProperty = dataSource.getContainerProperty(itemId, JobLogContainer.OUTPUT_TEXT);
		
		if (outputProperty != null) {
			output = (String) outputProperty.getValue();
		}
		
		showTextWindow("Job output", output);
	}

	public void showErrorOutput(Object itemId) {
		String error = "";
		Property<?> errorProperty = dataSource.getContainerProperty(itemId, JobLogContainer.ERROR_MESSAGE);
		
		if (errorProperty != null) {
			error = (String) errorProperty.getValue();
		}

		showTextWindow("Error message", error);
	}
	
	private void showTextWindow(String caption, String content) {
		
		Label textComponent = new Label(content);
		textComponent.setContentMode(ContentMode.PREFORMATTED);
		
		Window subWindow = new Window(caption);
		subWindow.setContent(textComponent);
		
		subWindow.setWidth(70, Unit.PERCENTAGE);
		subWindow.setHeight(90, Unit.PERCENTAGE);
		subWindow.center();
		
		this.getUI().addWindow(subWindow);
	}
}
