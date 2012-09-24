package fi.csc.microarray.manager.web.hbncontainer;

import java.util.Date;
import java.util.Random;

import org.hibernate.SessionFactory;
import org.hibernate.cfg.AnnotationConfiguration;
import org.hibernate.cfg.Environment;
import org.hibernate.classic.Session;
import org.hibernate.dialect.H2Dialect;

import fi.csc.microarray.manager.web.data.JobLogEntry;
import fi.csc.microarray.manager.web.util.RandomUtil;

public class JobLogHibernateUtil {

    private static final SessionFactory sessionFactory;
    //private static Type defaultType;

    static {
        try {
            AnnotationConfiguration cnf = new AnnotationConfiguration();
            
            //CREATE TABLE IF NOT EXISTS jobs (id VARCHAR(100) PRIMARY KEY, operation VARCHAR(200), status VARCHAR(200), starttime DATETIME DEFAULT NULL, endtime DATETIME DEFAULT NULL, 	wallclockTime INT DEFAULT NULL, errorMessage TEXT DEFAULT NULL, outputText TEXT DEFAULT NULL, username VARCHAR(200), compHost VARCHAR(500));JOBS ID JOBS
            //insert into jobs values ('8f99df7a-5c38-439f-8f7b-15f8bed88eda', 'importseq.sadl', 'COMPLETED', '2010-07-08 13:40:39.0', '2010-07-08 JOBS 13:40:39.0', 0, null, null, 'demo', 'hippu2.csc.fi');
            
            cnf.setProperty(Environment.DRIVER, "org.h2.Driver");
            //cnf.setProperty(Environment.URL, "jdbc:h2:~/test");
            cnf.setProperty(Environment.URL, "jdbc:h2:tcp://localhost/test");
            cnf.setProperty(Environment.USER, "sa");
            cnf.setProperty(Environment.DIALECT, H2Dialect.class.getName());          
            cnf.setProperty(Environment.SHOW_SQL, "true");
            cnf.setProperty(Environment.CURRENT_SESSION_CONTEXT_CLASS, "thread");
            
//Useful for testing, but HAS TO BE DISABLED FOR PRODUCTION, otherwise data table may be dropped
//			cnf.setProperty(Environment.HBM2DDL_AUTO, "create");
            
//    		String dbDriver = configuration.getString("manager", "jdbc-driver");
//    		String dbUrl = configuration.getString("manager", "database-url");
//    		boolean startWebConsole = configuration.getBoolean("manager", "start-web-console");
//    		String dbUsername = configuration.getString("manager", "database-username");
//    	    String dbPassword = configuration.getString("manager", "database-password");
//    	    int webConsolePort = configuration.getInt("manager", "web-console-port");
            
            cnf.addAnnotatedClass(JobLogEntry.class);

            sessionFactory = cnf.buildSessionFactory();

        } catch (Throwable ex) {
            // Make sure you log the exception, as it might be swallowed
            System.err.println("Initial SessionFactory creation failed." + ex);
            throw new ExceptionInInitializerError(ex);
        }
    }

    public static SessionFactory getSessionFactory() {
        return sessionFactory;
    }
    


    public static void insertExampleData(int count) {
        Session sess = getSessionFactory().getCurrentSession();
        sess.beginTransaction();

    	String[] states = new String[] { "COMPLETED", "FAILED" };
    	Date startTime;
    	Date endTime;
    	int wallclockTime;
//    	String errorMessage;
//    	String outputText;
    	
        Random rnd = new Random();

        JobLogEntry job;

        for (int i = 0; i < count; i++) {
            job = new JobLogEntry();
            
            job.setId("8f99dsdf7a-5c38-439f-8f7b-15f8bed" + Math.abs(rnd.nextLong()));
            job.setOperation(RandomUtil.getRandomOperation(rnd));
            
            Date start = RandomUtil.getRandomDate(rnd, 2009);
            job.setStartTime(start);
            int duration = rnd.nextInt(rnd.nextInt(600) + 1);
            Date end = (Date) start.clone();
            end.setMinutes(start.getMinutes() + duration);
            job.setEndTime(end);
            job.setWallclockTime(duration);
            
            if (rnd.nextDouble() < 0.02) {
            	job.setStatus(states[1]);
            	job.setErrorMessage(RandomUtil.getRandomError(rnd));
            } else {
            	job.setStatus(states[0]);
            }
            
            job.setOutputText(RandomUtil.getRandomOutputtext(rnd));
            
            job.setUsername(RandomUtil.getRandomUserName(rnd));
            job.setCompHost(RandomUtil.getRandomComp(rnd));

            sess.save(job);
        }
    }
}
