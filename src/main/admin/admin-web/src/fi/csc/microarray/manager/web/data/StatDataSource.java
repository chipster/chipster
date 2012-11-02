package fi.csc.microarray.manager.web.data;

import java.util.Calendar;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.hibernate.Session;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.transform.Transformers;

public class StatDataSource {
	
	public static final String ROW_COUNT = "count";

	
	public List<Map> getTopUsers(Session session) {
		
		//Start one year ago
		Calendar fromDate = new GregorianCalendar();
		fromDate.set(Calendar.YEAR, fromDate.get(Calendar.YEAR) - 1);
		
		@SuppressWarnings("unchecked")
		List<Map> results = session.createCriteria(JobLogEntry.class)
				.add(Restrictions.ge(JobLogContainer.START_TIME, fromDate.getTime()))
			    .setProjection( Projections.projectionList()
			        .add(Projections.rowCount(), ROW_COUNT)
			        .add(Projections.groupProperty(JobLogContainer.USERNAME), JobLogContainer.USERNAME)		        
			    )
			    .addOrder(Order.desc(ROW_COUNT))
			    .setMaxResults(10)
			    .setResultTransformer(Transformers.ALIAS_TO_ENTITY_MAP)
			    .list();
		
		//Create one empty row to carry the table column names, if the result was empty
		if (results == null || results.size() == 0) {	
			HashMap<String, Object> map = new HashMap<String, Object>();
			map.put(JobLogContainer.USERNAME, null);
			map.put(ROW_COUNT, null);
			results = new LinkedList<Map>();
			results.add(map);
		}
		
		return results;
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
		
		return results;
	}
	
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public List<Map> getJobsByMonth(Session session) {
		
		List results =  session.createQuery("select new map(year(startTime) as year, month(startTime) as month, count(*) as count) " +
				"from JobLogEntry " + 
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
		
		return results;
	}
	
	public Object[] getJobsByMonthColumnOrder() {
		return new Object[] { "year", "month", "count" };
	}

	@SuppressWarnings("rawtypes")
	public List<Map> getToolFails(Session session) {
		//Start one year ago
		Calendar fromDate = new GregorianCalendar();
		fromDate.set(Calendar.YEAR, fromDate.get(Calendar.YEAR) - 1);
		
		@SuppressWarnings("unchecked")
		List<Map> results = session.createCriteria(JobLogEntry.class)
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
		
		return results;
	}
	
	public Object[] getToolFailsColumnOrder() {
		return new Object[] { JobLogContainer.OPERATION, ROW_COUNT };
	}
}
