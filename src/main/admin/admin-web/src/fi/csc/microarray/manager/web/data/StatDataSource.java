package fi.csc.microarray.manager.web.data;

import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.hibernate.Criteria;
import org.hibernate.Session;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.ProjectionList;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.transform.Transformers;

public class StatDataSource {

	public static final String ROW_COUNT = "count";

	private TestAccountFilter testAccountFilter = new TestAccountFilter();

	@SuppressWarnings("rawtypes")
	public List<Map<Object, Object>> getTopUsers(Session session, boolean ignoreTestAccounts) {

		//Start one year ago
		Calendar fromDate = new GregorianCalendar();

		fromDate.set(Calendar.YEAR, fromDate.get(Calendar.YEAR) - 1);

		Criteria criteria = session.createCriteria(JobLogEntry.class);		
		testAccountFilter.addCriteriaForTestAccounts(session, ignoreTestAccounts, criteria);

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

	@SuppressWarnings("rawtypes")
	private Map getJobDateRange(Session session, boolean ignoreTestAccounts) {
		//Get min and max values of dates
		Criteria rangeCriteria = session.createCriteria(JobLogEntry.class);
		testAccountFilter.addCriteriaForTestAccounts(session, ignoreTestAccounts, rangeCriteria);
		ProjectionList projections = Projections.projectionList();
		projections.add(Projections.min(JobLogContainer.START_TIME), "min");
		projections.add(Projections.max(JobLogContainer.START_TIME), "max");
		rangeCriteria.setProjection(projections);
		rangeCriteria.setResultTransformer(Transformers.ALIAS_TO_ENTITY_MAP);

		@SuppressWarnings("unchecked")
		List<Map> rangeList = rangeCriteria.list();

		return rangeList.get(0);
	}

	public List<Map<Object, Object>> getMonthlyStats(Session session, boolean ignoreTestAccounts) {


		@SuppressWarnings("rawtypes")
		Map rangeMap = getJobDateRange(session, ignoreTestAccounts);

		Calendar minDate = new GregorianCalendar();
		Calendar maxDate = new GregorianCalendar();
		minDate.setTime((Date) rangeMap.get("min"));
		maxDate.setTime((Date) rangeMap.get("max"));

		int minYear = minDate.get(Calendar.YEAR);
		int maxYear = maxDate.get(Calendar.YEAR);
		int minMonth = minDate.get(Calendar.MONTH) + 1;
		int maxMonth = maxDate.get(Calendar.MONTH) + 1;

		List<Map<Object, Object>> results = new LinkedList<Map<Object, Object>>();

		for (int year = minYear; year <= maxYear; year++) {

			int month;

			if (year == minYear) {
				month = minMonth;
			} else {
				month = 1;
			}

			while (month <= 12) {

				if (year == maxYear && month > maxMonth) {
					break;
				}

				Calendar startDate = new GregorianCalendar();
				startDate.set(year, month - 1, 1, 0, 0, 0);

				Calendar endDate = (Calendar) startDate.clone();        		
				endDate.add(Calendar.MONTH, 1);

				Criteria criteria = session.createCriteria(JobLogEntry.class);		
				testAccountFilter.addCriteriaForTestAccounts(session, ignoreTestAccounts, criteria);

				criteria.add(Restrictions.ge(JobLogContainer.START_TIME, startDate.getTime()));
				criteria.add(Restrictions.lt(JobLogContainer.START_TIME, endDate.getTime()));

				criteria.setProjection( Projections.projectionList()
						.add(Projections.rowCount(), "jobCount")
						.add(Projections.countDistinct(JobLogContainer.USERNAME), "uniqueUsers")
						);
				criteria.setResultTransformer(Transformers.ALIAS_TO_ENTITY_MAP);        	

 		
				@SuppressWarnings({ "rawtypes", "unchecked" })
				List<Map> resultList = criteria.list();
				
				@SuppressWarnings("unchecked")
				Map<Object, Object> resultMap = resultList.get(0);

				resultMap.put("year", year);
				resultMap.put("month", month);

				results.add(resultMap);

				month++;
			}
		}

		return results;
	}

	public Object[] getMonthlyStatsColumnOrder() {
		return new Object[] { "year", "month", "uniqueUsers", "jobCount"};
	}

	public List<Map<Object, Object>> getYearlyStats(Session session, boolean ignoreTestAccounts) {

		@SuppressWarnings("rawtypes")
		Map rangeMap = getJobDateRange(session, ignoreTestAccounts);

		Calendar minDate = new GregorianCalendar();
		Calendar maxDate = new GregorianCalendar();
		minDate.setTime((Date) rangeMap.get("min"));
		maxDate.setTime((Date) rangeMap.get("max"));

		int minYear = minDate.get(Calendar.YEAR);
		int maxYear = maxDate.get(Calendar.YEAR);

		List<Map<Object, Object>> results = new LinkedList<Map<Object, Object>>();

		for (int year = minYear; year <= maxYear; year++) {

			Calendar startDate = new GregorianCalendar();
			startDate.set(year, 0, 1, 0, 0, 0);

			Calendar endDate = (Calendar) startDate.clone();        		
			endDate.add(Calendar.YEAR, 1);

			Criteria criteria = session.createCriteria(JobLogEntry.class);		
			testAccountFilter.addCriteriaForTestAccounts(session, ignoreTestAccounts, criteria);

			criteria.add(Restrictions.ge(JobLogContainer.START_TIME, startDate.getTime()));
			criteria.add(Restrictions.lt(JobLogContainer.START_TIME, endDate.getTime()));

			criteria.setProjection( Projections.projectionList()
					.add(Projections.rowCount(), "jobCount")
					.add(Projections.countDistinct(JobLogContainer.USERNAME), "uniqueUsers")
					);
			criteria.setResultTransformer(Transformers.ALIAS_TO_ENTITY_MAP);        	
 		
			@SuppressWarnings({ "rawtypes", "unchecked" })
			List<Map> resultList = criteria.list();

			@SuppressWarnings("unchecked")
			Map<Object, Object> resultMap = resultList.get(0);				
			resultMap.put("year", year);
			results.add(resultMap);
		}
		return results;
	}

	public Object[] getYearlyStatsColumnOrder() {
		return new Object[] { "year", "uniqueUsers", "jobCount"};
	}

	@SuppressWarnings("rawtypes")
	public List<Map<Object, Object>> getToolFails(Session session, boolean ignoreTestAccounts) {
		//Start one year ago
		Calendar fromDate = new GregorianCalendar();
		fromDate.set(Calendar.YEAR, fromDate.get(Calendar.YEAR) - 1);

		Criteria criteria = session.createCriteria(JobLogEntry.class);
		testAccountFilter.addCriteriaForTestAccounts(session, ignoreTestAccounts, criteria);

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
	
	@SuppressWarnings("rawtypes")
	public List<Map<Object, Object>> getToolUsage(Session session, boolean ignoreTestAccounts) {
		//Start one year ago
		Calendar fromDate = new GregorianCalendar();
		fromDate.set(Calendar.YEAR, fromDate.get(Calendar.YEAR) - 1);

		Criteria criteria = session.createCriteria(JobLogEntry.class);
		testAccountFilter.addCriteriaForTestAccounts(session, ignoreTestAccounts, criteria);

		@SuppressWarnings({"unchecked" })
		List<Map> results = criteria
		.add(Restrictions.ge(JobLogContainer.START_TIME, fromDate.getTime()))
		.setProjection( Projections.projectionList()
				.add(Projections.rowCount(), ROW_COUNT)
				.add(Projections.groupProperty(JobLogContainer.OPERATION), JobLogContainer.OPERATION)		        
				)
				.addOrder(Order.desc(ROW_COUNT))
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

	public Object[] getToolUsageColumnOrder() {
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
