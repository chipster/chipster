package fi.csc.microarray.manager.web.hbncontainer;

import org.hibernate.Session;

import com.vaadin.Application;
import com.vaadin.data.hbnutil.HbnContainer.SessionManager;
import com.vaadin.service.ApplicationContext.TransactionListener;

public class JobLogSessionManager implements SessionManager {
	
		private Application app;

		public JobLogSessionManager(Application app) {
			this.app = app;
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
	}