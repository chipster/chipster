package fi.csc.chipster.web.adminweb.ui;

import java.util.LinkedList;

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
import fi.csc.chipster.web.adminweb.data.JobLogContainer;

public class JobLogView extends VerticalLayout implements ClickListener, ValueChangeListener  {
	
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
		init();
	}

	public void init() {
		
		//dataSourceWrapper has to be initialized here because of the transactionListener, so lets init everything else 
		//here as well (and not in the constructor like elsewhere)		
		
		// do this before data source is attached to avoid one data update
		setSizeFull();
		
		try {
			dataSource = new JobLogContainer(this);			
//			dataSource.init();
			
			table = new JobLogTable(this);		
			
			table.setContainerDataSource(dataSource);					
			
		} catch (GenericJDBCException e) {
			//FIXME Show exception message and hide or disable all database based content 
			return;
		}
		
		table.setVisibleColumns(JobLogContainer.NATURAL_COL_ORDER);
		table.setColumnHeaders(JobLogContainer.COL_HEADERS_ENGLISH);
		
		table.setSortAscending(false);
		table.setSortContainerPropertyId(JobLogContainer.START_TIME);
		
		this.addComponent(getToolbar());
		this.addComponent(table);

		this.setExpandRatio(table, 1);
	}


	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			
			toolbarLayout = new HorizontalLayout();
			
			filterLayout = new HorizontalLayout();
			addFilter();
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

					applyFilters();
				}
			});
			
			/*
			 * Don't filter test accounts by default. Current way of filtering
			 * strings is slow, because H2 doesn't use index for these SQL
			 * queries (WHERE NOT username = '').
			 */
			ignoreTestAccounts = new CheckBox("Ignore test accounts", false);
			ignoreTestAccounts.addStyleName("toolbar-component");
			toolbarLayout.addComponent(ignoreTestAccounts);
			dataSource.setIgnoreTestAccounts(ignoreTestAccounts.getValue());
						
			ignoreTestAccounts.addValueChangeListener(new ValueChangeListener() {

				@Override
				public void valueChange(ValueChangeEvent arg0) {
					applyFilters();
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

	private void addFilter() {
		if (filters == null) {
			filters = new LinkedList<JobLogFilter>();
		}
		JobLogFilter filter = new JobLogFilter(this);
		filters.add(filter);
		filterLayout.addComponent(filter);
	}

	public void buttonClick(ClickEvent event) {
		final Button source = event.getButton();

		if (source == addFilterButton) {
			addFilter();
		}
	}

	public void valueChange(ValueChangeEvent event) {
		Property<?> property = event.getProperty();
		if (property == table) {
			//			Item item = personList.getItem(personList.getValue());
			//			if (item != personForm.getItemDataSource()) {
			//				personForm.setItemDataSource(item);
			//			}
		}
	}

	public void applyFilters() {

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

	public void clearFilters(JobLogFilter filter) {

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
