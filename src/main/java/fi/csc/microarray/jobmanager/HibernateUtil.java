package fi.csc.microarray.jobmanager;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.hibernate.SessionFactory;
import org.hibernate.boot.registry.StandardServiceRegistryBuilder;
import org.hibernate.cfg.Environment;

import fi.csc.microarray.config.Configuration;

public class HibernateUtil {

	private static Logger logger = Logger.getLogger(HibernateUtil.class.getName());
	
    private SessionFactory sessionFactory;

    public void buildSessionFactory(List<Class<?>> hibernateClasses, Configuration config) {
    	
    	
    	try {    		
    		
    		final org.hibernate.cfg.Configuration hibernateConf = new org.hibernate.cfg.Configuration();

    		String dbDriver = config.getString("jobmanager", "hibernate-driver");
    		String dbDialect = config.getString("jobmanager", "hibernate-dialect");
    		String dbUrl = config.getString("jobmanager", "hibernate-url");
    		String dbUsername = config.getString("jobmanager", "hibernate-username");
    		String dbPassword = config.getString("jobmanager", "hibernate-password");		
    		String showSql = config.getBoolean("jobmanager", "hibernate-show-sql") ? "true" : "false";
    		String dbSchemaUpdate = config.getString("jobmanager", "hibernate-schema");

    		hibernateConf.setProperty(Environment.DRIVER, dbDriver);
    		hibernateConf.setProperty(Environment.URL, dbUrl);
    		hibernateConf.setProperty(Environment.USER, dbUsername);
    		hibernateConf.setProperty(Environment.PASS, dbPassword);
    		hibernateConf.setProperty(Environment.DIALECT, dbDialect);
    		hibernateConf.setProperty(Environment.SHOW_SQL, showSql);
    		hibernateConf.setProperty(Environment.CURRENT_SESSION_CONTEXT_CLASS, "thread");
    		hibernateConf.setProperty("hibernate.hbm2ddl.auto", dbSchemaUpdate);    		
    		hibernateConf.setProperty("hibernate.c3p0.min_size", "3");
    		
    		for (Class<?> c : hibernateClasses) {
    			hibernateConf.addAnnotatedClass(c);
    		}    		    	   
    		
    		sessionFactory = hibernateConf.buildSessionFactory(
    				new StandardServiceRegistryBuilder()
    				.applySettings(hibernateConf.getProperties())
    				.build());
 
    	} catch (Throwable ex) {
    		logger.log(Level.SEVERE, "sessionFactory creation failed.", ex);
    		throw new ExceptionInInitializerError(ex);
    	}
    }

    public SessionFactory getSessionFactory() {
        return sessionFactory;
    }

	public org.hibernate.Session beginTransaction() {
		org.hibernate.Session session = getSessionFactory().getCurrentSession();
		session.beginTransaction();
		return session;
	}

	public void commit() {
		getSessionFactory().getCurrentSession().getTransaction().commit();
	}
	
	public void rollback() {
		getSessionFactory().getCurrentSession().getTransaction().rollback();
	}

	public org.hibernate.Session session() {
		return getSessionFactory().getCurrentSession();
	}

	public void rollbackIfActive() {
		if (session().getTransaction().isActive()) {
			session().getTransaction().rollback();
		}		
	}
}