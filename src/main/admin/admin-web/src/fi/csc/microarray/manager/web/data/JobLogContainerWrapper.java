package fi.csc.microarray.manager.web.data;

import org.hibernate.Session;

import com.vaadin.Application;
import com.vaadin.data.Container;
import com.vaadin.data.Property;
import com.vaadin.data.hbnutil.ContainerFilter;
import com.vaadin.data.hbnutil.HbnContainer;
import com.vaadin.data.hbnutil.HbnContainer.SessionManager;
import com.vaadin.service.ApplicationContext.TransactionListener;

import fi.csc.microarray.manager.web.hbncontainer.JobLogHibernateUtil;
import fi.csc.microarray.manager.web.ui.JobLogView;


public class JobLogContainerWrapper implements SessionManager {
	
	/**
	 * Natural property order for Service bean. Used in tables and forms.
	 */
	public static final Object[] NATURAL_COL_ORDER = new Object[] {
		"username", "operation", "status", "compHost", "startTime", "endTime", "wallclockTime", "outputLink", "errorLink"}; //, "errorMessage", "outputText" };

	/**
	 * "Human readable" captions for properties in same order as in
	 * NATURAL_COL_ORDER.
	 */
	public static final String[] COL_HEADERS_ENGLISH = new String[] {
		"Username", "Operation", "Status", "Comp host", "Start time", "End time", "Wall clock time", "", "Error"}; //, "Error message", "Output text" };



	private HbnContainer<JobLogEntry> hbnContainer;
	private Application app;
	
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
		Session sess = JobLogHibernateUtil.getSessionFactory().getCurrentSession();
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
		Session currentSession = JobLogHibernateUtil.getSessionFactory()
				.getCurrentSession();
		if (!currentSession.getTransaction().isActive()) {
			currentSession.beginTransaction();
		}
		return currentSession;
	}


	public JobLogContainerWrapper(JobLogView view) {
		
		this.app = view.getApp();
	}

	public void init() {
		hbnContainer = new HbnContainer<JobLogEntry>(JobLogEntry.class, this);
		
		attachVaadinTransactionListener();
		
	}

	public Container getDataSource() {
		return hbnContainer;
	}

	public void removeAllContainerFilters() {
		hbnContainer.removeAllContainerFilters();
	}

	public void addContainerFilter(ContainerFilter containerFilter) {
		hbnContainer.addContainerFilter(containerFilter);
	}

	public Property getContainerProperty(Object itemId, String propertyId) {
		return hbnContainer.getContainerProperty(itemId, propertyId);
	}

}
