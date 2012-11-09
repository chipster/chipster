package fi.csc.microarray.manager.web.data;

import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.hibernate.Session;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.transform.Transformers;
import org.springframework.web.util.HtmlUtils;

public class StatDataSource {
	
	public static final String ROW_COUNT = "count";

	
	@SuppressWarnings("rawtypes")
	public List<Map<Object, Object>> getTopUsers(Session session) {
		
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
		
		return htmlEscape(results);
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
	public List<Map<Object, Object>> getJobsByMonth(Session session) {
		
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
		
		return htmlEscape(results);
	}
	
	public Object[] getJobsByMonthColumnOrder() {
		return new Object[] { "year", "month", "count" };
	}

	@SuppressWarnings("rawtypes")
	public List<Map<Object, Object>> getToolFails(Session session) {
		//Start one year ago
		Calendar fromDate = new GregorianCalendar();
		fromDate.set(Calendar.YEAR, fromDate.get(Calendar.YEAR) - 1);
		
		@SuppressWarnings({"unchecked" })
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
		
		return htmlEscape(results);
	}
	
	public Object[] getToolFailsColumnOrder() {
		return new Object[] { JobLogContainer.OPERATION, ROW_COUNT };
	}
		
	private List<Map<Object, Object>> htmlEscape(@SuppressWarnings("rawtypes") List<Map> unescapedMaps) {
		
		List<Map<Object, Object>> escapedMaps = new LinkedList<Map<Object, Object>>();
		
		for (@SuppressWarnings("rawtypes") Map unescapedMap : unescapedMaps) {
			
			Map<Object, Object> escapedMap = new HashMap<Object, Object>();
			
			for (Object entryObj  : unescapedMap.entrySet()) {
				
				@SuppressWarnings("rawtypes")
				Entry entry = (Entry)entryObj;
				
				Object key = entry.getKey();
				Object value = entry.getValue();
				
				if (key instanceof String) {
					key = HtmlUtils.htmlEscape((String)key);
				}
				
				if (value instanceof String) {
					value = HtmlUtils.htmlEscape((String)value);
				}
				
				escapedMap.put(key, value);
			}
			escapedMaps.add(escapedMap);
		}
		
		return escapedMaps;
	}
}
