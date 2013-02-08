package fi.csc.microarray.manager.web.data;

import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.hibernate.Criteria;
import org.hibernate.Session;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.transform.Transformers;

public class StatDataSource {

	public static final String ROW_COUNT = "count";

	private List<String> testAccounts;

	public List<String> getTestAccounts(Session session) {

		if (testAccounts == null) {
			@SuppressWarnings("unchecked")
			List<AccountEntry> results = session.createCriteria(AccountEntry.class)
			//FIXME
			//.add(Restrictions.eq("ignoreInStatistics", true))
			.list();

			testAccounts = new LinkedList<String>();

			for (AccountEntry entry : results) {
				testAccounts.add(entry.getUsername());
			}
		}
		return testAccounts;
	}
	
	private void addCriteriaForTestAccounts(Session session,
			boolean ignoreTestAccounts, Criteria criteria) {
		
		if (ignoreTestAccounts) {
			for (String account : getTestAccounts(session)) {
				criteria.add(Restrictions.not(Restrictions.eq(JobLogContainer.USERNAME, account)));
			}
		}
	}
	
	private String getHqlForTestAccounts(Session session, boolean ignoreTestAccounts) {
		
		String hql = "where not ";
		
		if (ignoreTestAccounts) {
			
			boolean isFirst = true;
			
			List<String> accounts = getTestAccounts(session);			
			
			for (String account : accounts) {
				
				//Check that we don't add anything dangerous to the query (only alphanumeric and not empty)
				if (account.matches("^[a-zA-Z0-9]+$")) {
					
					if (!isFirst) {
						hql += "and not ";
					}
					isFirst = false;
					
					hql += JobLogContainer.USERNAME + "='" + account + "' ";
				}
			}
			
			if (accounts.size() > 0) {
				
				return hql;
			}
		}
		
		return "";
	}


	@SuppressWarnings("rawtypes")
	public List<Map<Object, Object>> getTopUsers(Session session, boolean ignoreTestAccounts) {

		//Start one year ago
		Calendar fromDate = new GregorianCalendar();
		
		fromDate.set(Calendar.YEAR, fromDate.get(Calendar.YEAR) - 1);

		Criteria criteria = session.createCriteria(JobLogEntry.class);		
		addCriteriaForTestAccounts(session, ignoreTestAccounts, criteria);
		
		criteria.add(Restrictions.ge(JobLogContainer.START_TIME, fromDate.getTime()));
		criteria.setProjection( Projections.projectionList()
				.add(Projections.rowCount(), ROW_COUNT)
				.add(Projections.groupProperty(JobLogContainer.USERNAME), JobLogContainer.USERNAME)		        
		);
		criteria.addOrder(Order.desc(ROW_COUNT));
		criteria.setMaxResults(10);
		criteria.setResultTransformer(Transformers.ALIAS_TO_ENTITY_MAP);
		
		@SuppressWarnings("unchecked") 		
		List<Map> results = criteria.list();

		//Create one empty row to carry the table column names, if the result was empty
		if (results == null || results.size() == 0) {	
			HashMap<String, Object> map = new HashMap<String, Object>();
			map.put(JobLogContainer.USERNAME, null);
			map.put(ROW_COUNT, null);
			results = new LinkedList<Map>();
			results.add(map);
		}

		return typedMap(results);
	}

	public Object[] getTopUsersColumnOrder() {
		return new Object[] { JobLogContainer.USERNAME, ROW_COUNT };
	}

	public List<JobLogEntry> getLatestJobs(Session session) {

		@SuppressWarnings("unchecked")
		List<JobLogEntry> results = session.createCriteria(JobLogEntry.class)
		.addOrder(Order.desc(JobLogContainer.START_TIME))
		.setMaxResults(10)
		.list();

		//JobLogEntry does htmlEscaping
		return results;
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	public List<Map<Object, Object>> getJobsByMonth(Session session, boolean ignoreTestAccounts) {

		List results =  session.createQuery("select new map(year(startTime) as year, month(startTime) as month, count(*) as count) " +
				"from JobLogEntry " + 
				getHqlForTestAccounts(session, ignoreTestAccounts)  +
				"group by year(startTime), month(startTime) " +
				"order by year(startTime) asc, month(startTime) asc").list();		

		//Create one empty row to carry the table column names, if the result was empty
		if (results == null || results.size() == 0) {	
			HashMap<String, Object> map = new HashMap<String, Object>();
			Object[] cols = getJobsByMonthColumnOrder();
			map.put((String) cols[0], null);
			map.put((String) cols[1], null);
			map.put((String) cols[2], null);
			results = new LinkedList<Map>();
			results.add(map);
		}

		return typedMap(results);
	}

	public Object[] getJobsByMonthColumnOrder() {
		return new Object[] { "year", "month", "count" };
	}

	@SuppressWarnings("rawtypes")
	public List<Map<Object, Object>> getToolFails(Session session, boolean ignoreTestAccounts) {
		//Start one year ago
		Calendar fromDate = new GregorianCalendar();
		fromDate.set(Calendar.YEAR, fromDate.get(Calendar.YEAR) - 1);

		Criteria criteria = session.createCriteria(JobLogEntry.class);
		addCriteriaForTestAccounts(session, ignoreTestAccounts, criteria);
		
		@SuppressWarnings({"unchecked" })
		List<Map> results = criteria
		.add(Restrictions.ge(JobLogContainer.START_TIME, fromDate.getTime()))
		.add(Restrictions.eq(JobLogContainer.STATUS, JobLogContainer.STATUS_FAIL_VALUE))
		.setProjection( Projections.projectionList()
				.add(Projections.rowCount(), ROW_COUNT)
				.add(Projections.groupProperty(JobLogContainer.OPERATION), JobLogContainer.OPERATION)		        
				)
				.addOrder(Order.desc(ROW_COUNT))
				.setMaxResults(10)
				.setResultTransformer(Transformers.ALIAS_TO_ENTITY_MAP)
				.list();

		//Create one empty row to carry the table column names, if the result was empty
		if (results == null || results.size() == 0) {	
			HashMap<String, Object> map = new HashMap<String, Object>();
			map.put(JobLogContainer.OPERATION, null);
			map.put(ROW_COUNT, null);
			results = new LinkedList<Map>();
			results.add(map);
		}

		return typedMap(results);
	}

	public Object[] getToolFailsColumnOrder() {
		return new Object[] { JobLogContainer.OPERATION, ROW_COUNT };
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	private List<Map<Object, Object>> typedMap(List<Map> rawtype) {

		List<Map<Object, Object>> typed = new LinkedList<Map<Object, Object>>();

		for (Map rawItem : rawtype) {			
			typed.add(rawItem);
		}

		return typed;
	}
}
