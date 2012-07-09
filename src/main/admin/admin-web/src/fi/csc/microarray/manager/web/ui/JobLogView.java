package fi.csc.microarray.manager.web.ui;

import java.util.LinkedList;

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
import com.vaadin.ui.Window;
import com.vaadin.ui.Window.Notification;

import fi.csc.microarray.manager.web.ChipsterAdminApplication;
import fi.csc.microarray.manager.web.data.JobLogContainer;
import fi.csc.microarray.manager.web.hbncontainer.JobLogHibernateUtil;
import fi.csc.microarray.manager.web.hbncontainer.JobLogSessionManager;

public class JobLogView extends VerticalLayout implements ClickListener, ValueChangeListener  {
	
	private HorizontalLayout toolbarLayout;

	private Button refreshButton = new Button("Insert 1k rows");
	private Button addSearchButton = new Button();

	private JobLogTable table;
	private JobLogContainer dataSource;

	private ChipsterAdminApplication app;

	private LinkedList<JobLogSearch> searches;


	private HorizontalLayout searchLayout;


	public JobLogView(ChipsterAdminApplication app) {

		this.app = app;
	}

	public void init() {
		
		//dataSourceWrapper has to be initialized here because of the transactionListener, so lets init everything else 
		//here as well (and not in the constructor like elsewhere)
		
		dataSource = new JobLogContainer(this, new JobLogSessionManager(app)); 
		dataSource.init();
				
		table = new JobLogTable(this);
		table.setContainerDataSource(dataSource);

		table.setVisibleColumns(JobLogContainer.NATURAL_COL_ORDER);
		table.setColumnHeaders(JobLogContainer.COL_HEADERS_ENGLISH);
		
		this.addComponent(getToolbar());
		this.addComponent(table);

		setSizeFull();
		this.setExpandRatio(table, 1);
	}


	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			
			toolbarLayout = new HorizontalLayout();
			

			
			refreshButton.addListener((ClickListener)this);
			refreshButton.setIcon(new ThemeResource("../runo/icons/32/document-add.png"));
			toolbarLayout.addComponent(refreshButton);
			
			searchLayout = new HorizontalLayout();
			addSearch();
			toolbarLayout.addComponent(searchLayout);
			addSearchButton.addListener((ClickListener)this);
			addSearchButton.setIcon(new ThemeResource("crystal/edit_add.png"));
			addSearchButton.setDescription("Add another search");
			addSearchButton.addStyleName("search-button");
			
			HorizontalLayout buttonBorder = new HorizontalLayout();
			buttonBorder.addStyleName("search-filter");
			buttonBorder.addComponent(addSearchButton);
			toolbarLayout.addComponent(buttonBorder);
			
			Button searchButton = new Button();
			searchButton.setIcon(new ThemeResource("crystal/mail_find.png"));
			searchButton.setDescription("Search");
			//searchButton.addStyleName("search-button");
			toolbarLayout.addComponent(searchButton);
			
			searchButton.addListener(new Button.ClickListener() {

				public void buttonClick(ClickEvent event) {

					performSearch();
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

	private void addSearch() {
		if (searches == null) {
			searches = new LinkedList<JobLogSearch>();
		}
		JobLogSearch search = new JobLogSearch(this);
		searches.add(search);
		searchLayout.addComponent(search);
	}

	public void buttonClick(ClickEvent event) {
		final Button source = event.getButton();

		if (source == refreshButton) {
			JobLogHibernateUtil.insertExampleData(1000);
		} else if (source == addSearchButton) {
			addSearch();
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

	public void performSearch() {

		updateContainerFilters();

		app.getMainWindow().showNotification(
				"Found "
						+ table.getContainerDataSource().size() + " item(s)",
						Notification.TYPE_TRAY_NOTIFICATION);
	}
	
	private void updateContainerFilters() {
		dataSource.removeAllContainerFilters();

		for (JobLogSearch iteratedSearch : searches) {
			if (iteratedSearch.getContainerFilter() != null) {
				dataSource.addContainerFilter(iteratedSearch.getContainerFilter());
			}
		}
	}

	public void clearSearch(JobLogSearch search) {

		if (searches.size() > 1) {
			searchLayout.removeComponent(search);
			searches.remove(search);
		} else {
			searches.get(0).clear();
		}

		//It's not possible to remove just one filter, so let's do it in the hard way
		updateContainerFilters();
	}

	public void showOutput(Object itemId) {
		
		String output = "";
		Property outputProperty = dataSource.getContainerProperty(itemId, JobLogContainer.OUTPUT_TEXT);
		
		if (outputProperty != null) {
			output = (String) outputProperty.getValue();
		}
		
		showTextWindow("Job output", output);
	}

	public void showErrorOutput(Object itemId) {
		String error = "";
		Property errorProperty = dataSource.getContainerProperty(itemId, JobLogContainer.ERROR_MESSAGE);
		
		if (errorProperty != null) {
			error = (String) errorProperty.getValue();
		}

		showTextWindow("Error message", error);
	}
	
	private void showTextWindow(String caption, String content) {
		
		Label textComponent = new Label(content);
		textComponent.setContentMode(Label.CONTENT_PREFORMATTED);
		
		Window subWindow = new Window(caption);
		subWindow.addComponent(textComponent);
		
		subWindow.setWidth(70, UNITS_PERCENTAGE);
		subWindow.setHeight(90, UNITS_PERCENTAGE);
		subWindow.center();
		
		this.getWindow().addWindow(subWindow);
	}
}
