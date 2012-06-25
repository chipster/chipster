package fi.csc.microarray.manager.web.ui;

import java.util.LinkedList;

import org.hibernate.Session;

import com.vaadin.Application;
import com.vaadin.data.Property;
import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.data.hbnutil.HbnContainer;
import com.vaadin.data.hbnutil.HbnContainer.SessionManager;
import com.vaadin.service.ApplicationContext.TransactionListener;
import com.vaadin.terminal.ThemeResource;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.Table;
import com.vaadin.ui.VerticalLayout;
import com.vaadin.ui.Window.Notification;

import fi.csc.microarray.manager.web.ChipsterAdminApplication;
import fi.csc.microarray.manager.web.data.JobLogEntry;
import fi.csc.microarray.manager.web.hbncontainer.HibernateUtil;

public class JobLogView extends VerticalLayout implements ClickListener, ValueChangeListener, SessionManager {
	
	
	

	/**
	 * Natural property order for Service bean. Used in tables and forms.
	 */
	public static final Object[] NATURAL_COL_ORDER = new Object[] {
		"username", "operation", "state", "compHost", "startTime", "endTime", "wallclockTime", "errorMessage", "outputText" };

	/**
	 * "Human readable" captions for properties in same order as in
	 * NATURAL_COL_ORDER.
	 */
	public static final String[] COL_HEADERS_ENGLISH = new String[] {
		"Username", "Operation", "State", "Comp host", "Start time", "End time", "Wall clock time", "Error message", "Output text" };


	private HorizontalLayout toolbarLayout;

	private Button refreshButton = new Button("Insert 160k rows");
	private Button addSearchButton = new Button();

	private Table table;
	private HbnContainer<JobLogEntry> dataSource;

	private ChipsterAdminApplication app;

	private LinkedList<JobLogSearch> searches;


	private HorizontalLayout searchLayout;


	public JobLogView(ChipsterAdminApplication app) {

		this.app = app;

		buildView();
	}

	public void init() {
		attachVaadinTransactionListener();
	}


	/**
	 * HbnContainer: We are using session-per-request pattern with Hibernate. By using
	 * Vaadin's transaction listener we can easily ensure that session is closed
	 * on each request without polluting our program code with extra logic.
	 */
	public void attachVaadinTransactionListener() {
		app.getContext().addTransactionListener(new TransactionListener() {
			public void transactionEnd(Application application,
					Object transactionData) {
				// Transaction listener gets fired for all (Http) sessions
				// of Vaadin applications, checking to be this one.
				if (application == app) {
					closeSession();
				}
			}

			public void transactionStart(Application application,
					Object transactionData) {

			}
		});
	}

	/**
	 * HbnContainer
	 */
	private void closeSession() {
		Session sess = HibernateUtil.getSessionFactory().getCurrentSession();
		if (sess.getTransaction().isActive()) {
			sess.getTransaction().commit();
		}
		if (sess.isOpen()) {
			sess.close();
		}
	}

	/**
	 * HbnContainer: Used to get current Hibernate session. Also ensures an open Hibernate
	 * transaction.
	 */
	public Session getSession() {
		Session currentSession = HibernateUtil.getSessionFactory()
				.getCurrentSession();
		if (!currentSession.getTransaction().isActive()) {
			currentSession.beginTransaction();
		}
		return currentSession;
	}

	/**
	 * Builds a simple view for application with Table and a row of buttons
	 * below it.
	 */
	private void buildView() {

		populateAndConfigureTable();

		this.addComponent(getToolbar());
		this.addComponent(table);

		setSizeFull();
		table.setSizeFull();
		this.setExpandRatio(table, 1);
	}

	protected void populateAndConfigureTable() {
		table = new Table();

		table.setWidth("100%");
		table.setSelectable(true);
		table.setImmediate(true);

		loadJobLog();

		table.setVisibleColumns(NATURAL_COL_ORDER);
		table.setColumnHeaders(COL_HEADERS_ENGLISH);
	}

	/**
	 * Loads container to Table
	 */
	protected void loadJobLog() { 

		dataSource = new HbnContainer<JobLogEntry>(JobLogEntry.class, this);

		table.setContainerDataSource(dataSource);
	}


	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			
			toolbarLayout = new HorizontalLayout();
			
			Label spaceEater = new Label(" ");
			toolbarLayout.addComponent(spaceEater);
			toolbarLayout.setExpandRatio(spaceEater, 1);
			
//			refreshButton.addListener((ClickListener)this);
//			refreshButton.setIcon(new ThemeResource("../runo/icons/32/document-add.png"));
//			toolbarLayout.addComponent(refreshButton);
			
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
			HibernateUtil.insertExampleData(160000);
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
}
