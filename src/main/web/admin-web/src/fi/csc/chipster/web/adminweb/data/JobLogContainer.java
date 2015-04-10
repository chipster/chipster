package fi.csc.chipster.web.adminweb.data;


import org.hibernate.Criteria;
import org.hibernate.Session;

import com.vaadin.data.hbnutil.HbnContainer;

import fi.csc.chipster.web.adminweb.hbncontainer.HibernateUtil;
import fi.csc.chipster.web.adminweb.ui.JobLogView;


public class JobLogContainer extends HbnContainer<JobLogEntry> {

	public static final String USERNAME = "username";
	public static final String OPERATION = "operation";
	public static final String STATUS = "status";
	public static final String COMPHOST = "compHost";
	public static final String START_TIME = "startTime";
	public static final String END_TIME = "endTime";
	public static final String WALLCLOCK_TIME = "wallclockTime";
	public static final String ERROR_MESSAGE = "errorMessage";
	public static final String OUTPUT_TEXT = "outputText";
	public static final String OUTPUT_LINK = "outputLink";
	public static final String ERROR_LINK = "errorLink";

	public static final Object[] NATURAL_COL_ORDER  = new String[] {
		USERNAME, 		OPERATION, 		COMPHOST, 		WALLCLOCK_TIME, 	END_TIME, 	OUTPUT_LINK, 	ERROR_LINK, 	STATUS	};

	public static final String[] COL_HEADERS_ENGLISH = new String[] {
		"Username", 	"Operation", 	"Comp host", 	"Wall clock time",	"End time",	"", 			"Error", 		"Status" };
	
	
	public static final String STATUS_FAIL_VALUE = "FAILED";
	
	private TestAccountFilter testAccountFilter = new TestAccountFilter();
	private boolean ignoreTestAccounts;

	public JobLogContainer(JobLogView view) {

		super(JobLogEntry.class, HibernateUtil.getSessionFactory());
	}
	
	@Override
	protected Criteria getBaseCriteria() {
		Session session = HibernateUtil.getSessionFactory().getCurrentSession();
		
		//Filter test accounts
		Criteria criteria = super.getBaseCriteria();		
		testAccountFilter.addCriteriaForTestAccounts(session, ignoreTestAccounts, criteria);
		
		return criteria;
	}

	public void setIgnoreTestAccounts(boolean ignoreTestAccounts) {
		this.ignoreTestAccounts = ignoreTestAccounts;
	}
}
