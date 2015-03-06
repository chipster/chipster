package fi.csc.chipster.web.adminweb.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import org.hibernate.Criteria;
import org.hibernate.Session;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.ProjectionList;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.transform.Transformers;
import org.springframework.core.io.ClassPathResource;

import com.vaadin.server.FileResource;
import com.vaadin.server.VaadinService;

public class StatDataSource {

	public static final String ROW_COUNT = "count";
	public static final String YEAR = "year";
	public static final String MONTH = "month";
	public static final String JOB_COUNT = "jobCount";
	public static final String UNIQUE_USERS = "uniqueUsers";
	public static final String MICROARRAY = "microarray";
	public static final String NGS = "ngs";
	public static final String UNRESOLVED = "unresolved";
	public static final String OTHER = "other";
	
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
		criteria.setMaxResults(1000);
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

		Date min = (Date) rangeMap.get("min");
		Date max = (Date) rangeMap.get("max");
		
		List<Map<Object, Object>> results = new LinkedList<Map<Object, Object>>();
		
		if (min != null || max != null) {
			
			Calendar minDate = new GregorianCalendar();
			Calendar maxDate = new GregorianCalendar();
			minDate.setTime(min);
			maxDate.setTime(max);


			int minYear = minDate.get(Calendar.YEAR);
			int maxYear = maxDate.get(Calendar.YEAR);
			int minMonth = minDate.get(Calendar.MONTH) + 1;
			int maxMonth = maxDate.get(Calendar.MONTH) + 1;


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
							.add(Projections.rowCount(), JOB_COUNT)
							.add(Projections.countDistinct(JobLogContainer.USERNAME), UNIQUE_USERS)
							);
					criteria.setResultTransformer(Transformers.ALIAS_TO_ENTITY_MAP);        	


					@SuppressWarnings({ "rawtypes", "unchecked" })
					List<Map> resultList = criteria.list();

					@SuppressWarnings("unchecked")
					Map<Object, Object> resultMap = resultList.get(0);

					resultMap.put(YEAR, year);
					resultMap.put(MONTH, month);

					results.add(resultMap);

					month++;
				}
			}
			
		} else {
			
			Map<Object, Object> resultMap = new HashMap<Object, Object>();
				
			resultMap.put(YEAR, null);
			resultMap.put(MONTH, null);
			resultMap.put(UNIQUE_USERS, null);
			resultMap.put(JOB_COUNT, null);

			results.add(resultMap);
		}
		
		return results;
	}

	public Object[] getMonthlyStatsColumnOrder() {
		return new Object[] { YEAR, MONTH, UNIQUE_USERS, JOB_COUNT};
	}

	public List<Map<Object, Object>> getYearlyStats(Session session, boolean ignoreTestAccounts) {

		@SuppressWarnings("rawtypes")
		Map rangeMap = getJobDateRange(session, ignoreTestAccounts);

		Calendar minDate = new GregorianCalendar();
		Calendar maxDate = new GregorianCalendar();
		
		Date min = (Date) rangeMap.get("min");
		Date max = (Date) rangeMap.get("max");
		
		List<Map<Object, Object>> results = new LinkedList<Map<Object, Object>>();
		
		if (min != null || max != null) {
			
			minDate.setTime(min);
			maxDate.setTime(max);

			int minYear = minDate.get(Calendar.YEAR);
			int maxYear = maxDate.get(Calendar.YEAR);

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
						.add(Projections.rowCount(), JOB_COUNT)
						.add(Projections.countDistinct(JobLogContainer.USERNAME), UNIQUE_USERS)
						);
				criteria.setResultTransformer(Transformers.ALIAS_TO_ENTITY_MAP);        	

				@SuppressWarnings({ "rawtypes", "unchecked" })
				List<Map> resultList = criteria.list();

				@SuppressWarnings("unchecked")
				Map<Object, Object> resultMap = resultList.get(0);				
				resultMap.put(YEAR, year);
				results.add(resultMap);
			}
		} else {
			
			Map<Object, Object> resultMap = new HashMap<Object, Object>();
			
			for (Object key : this.getYearlyStatsColumnOrder()) {	
				resultMap.put(key, null);
			}
			
			results.add(resultMap);
		}
		return results;
	}

	public Object[] getYearlyStatsColumnOrder() {
		return new Object[] { YEAR, UNIQUE_USERS, JOB_COUNT };
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
				//.setMaxResults(10)
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
	
	public List<Map<Object, Object>> getModuleUsage(Session session, boolean ignoreTestAccounts) {

		//Get date of first and last job
		@SuppressWarnings("rawtypes")
		Map rangeMap = getJobDateRange(session, ignoreTestAccounts);

		Calendar minDate = new GregorianCalendar();
		Calendar maxDate = new GregorianCalendar();
		
		Date min = (Date) rangeMap.get("min");
		Date max = (Date) rangeMap.get("max");
		
		List<Map<Object, Object>> moduleResults = new LinkedList<Map<Object, Object>>();	
		
		if (min != null || max != null) {
			
			minDate.setTime(min);
			maxDate.setTime(max);
		
			int minYear = minDate.get(Calendar.YEAR);
			int maxYear = maxDate.get(Calendar.YEAR);

			List<Map<Object, Object>> results = new LinkedList<Map<Object, Object>>();
			
			//Get a yearly job count for each tool 		
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
						.add(Projections.rowCount(), ROW_COUNT)
						.add(Projections.groupProperty(JobLogContainer.OPERATION), JobLogContainer.OPERATION)		        
						);

				criteria.setResultTransformer(Transformers.ALIAS_TO_ENTITY_MAP);        	

				@SuppressWarnings({ "rawtypes", "unchecked" })
				List<Map> resultList = criteria.list();

				for (@SuppressWarnings("rawtypes") Map item : resultList) {
					@SuppressWarnings("unchecked")
					Map<Object, Object> toolCount = item;
					toolCount.put(YEAR, year);
					results.add(toolCount);
				}
			}	

			String microarray = null;
			String ngs = null;

			try {
				//Try to load files from jar (in case of real server)
				microarray = readFile(new ClassPathResource("WebContent/WEB-INF/microarray-module-copy.xml").getInputStream());
				ngs = readFile(new ClassPathResource("WebContent/WEB-INF/ngs-module-copy.xml").getInputStream());

			} catch (FileNotFoundException e) {

				try {

					//This is not a real server, because files were not found from jar. Use real files instead (in case of running directly from sources)
					String basepath = VaadinService.getCurrent().getBaseDirectory().getAbsolutePath();

					microarray = readFile(new FileResource(new File(basepath + "/WEB-INF/microarray-module-copy.xml")).getSourceFile());
					ngs = readFile(new FileResource(new File(basepath + "/WEB-INF/ngs-module-copy.xml")).getSourceFile());

				} catch (IOException e1) {
					e1.printStackTrace();
				}

			} catch (IOException e1) {
				e1.printStackTrace();
			}				

			//Find out a module for each tool		
			Map<String, String> toolModuleMap = new HashMap<String, String>();

			for (Map<Object, Object> resultMap : results) {
				String tool = (String) resultMap.get(JobLogContainer.OPERATION);

				if (!toolModuleMap.containsKey(tool)) {

					boolean isMicroarray = microarray.contains(tool);
					boolean isNgs = ngs.contains(tool);

					if (isMicroarray && !isNgs) {
						toolModuleMap.put(tool, MICROARRAY);					
					} else if (isNgs && !isMicroarray) {
						toolModuleMap.put(tool, NGS);
					} else if (isNgs && isMicroarray) {
						toolModuleMap.put(tool, UNRESOLVED);
					} else if (!isNgs && !isMicroarray) {
						toolModuleMap.put(tool, OTHER);
					}
				}
			}				

			//Count yearly job count for each module
			Map<Integer, Map<String, Integer>> yearModuleCount = new TreeMap<Integer, Map<String, Integer>>();

			for (Map<Object, Object> resultMap : results) {
				String tool = (String) resultMap.get(JobLogContainer.OPERATION);
				Integer year = (Integer) resultMap.get(YEAR);
				Integer toolCount = (int)(long)(Long)resultMap.get(ROW_COUNT);

				if (!yearModuleCount.containsKey(year)) {
					yearModuleCount.put(year, new HashMap<String, Integer>());
				}

				Map<String, Integer> moduleCountMap = yearModuleCount.get(year);

				String module = toolModuleMap.get(tool);

				if (!moduleCountMap.containsKey(module)) {
					moduleCountMap.put(module, toolCount);
				} else {
					moduleCountMap.put(module, moduleCountMap.get(module) + toolCount);
				}				
			}		

			//Convert nested maps to list of maps	
			for (Entry<Integer, Map<String, Integer>> entry : yearModuleCount.entrySet()) {

				Map<Object, Object> moduleMap = new HashMap<Object, Object>();

				int microarrayCount = 0;
				int ngsCount = 0;
				int unresolvedCount = 0;
				int otherCount = 0;

				if (entry.getValue().containsKey(MICROARRAY)) {
					microarrayCount = entry.getValue().get(MICROARRAY);
				}
				if (entry.getValue().containsKey(NGS)) {
					ngsCount = entry.getValue().get(NGS);
				}
				if (entry.getValue().containsKey(UNRESOLVED)) {
					unresolvedCount = entry.getValue().get(UNRESOLVED);
				}
				if (entry.getValue().containsKey(OTHER)) {
					otherCount = entry.getValue().get(OTHER);
				}

				moduleMap.put(YEAR, entry.getKey());
				moduleMap.put(MICROARRAY, microarrayCount);
				moduleMap.put(NGS, ngsCount);
				moduleMap.put(UNRESOLVED, unresolvedCount);
				moduleMap.put(OTHER, otherCount);

				moduleResults.add(moduleMap);
			}
		} else {
			
			Map<Object, Object> moduleMap = new HashMap<Object, Object>();
			
			for (Object key : this.getModuleUsageColumnOrder()) {	
				moduleMap.put(key, null);
			}
			
			moduleResults.add(moduleMap);
		}
		return moduleResults;
	}

	public Object[] getModuleUsageColumnOrder() {
		return new Object[] { YEAR, MICROARRAY, NGS, UNRESOLVED, OTHER };
	}
	

	@SuppressWarnings({ "unchecked", "rawtypes" })
	private List<Map<Object, Object>> typedMap(List<Map> rawtype) {

		List<Map<Object, Object>> typed = new LinkedList<Map<Object, Object>>();

		for (Map rawItem : rawtype) {			
			typed.add(rawItem);
		}

		return typed;
	}

	private String readFile(File file) throws IOException {
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
		StringBuilder stringBuilder = new StringBuilder();
		String ls = System.getProperty("line.separator");

		while ((line = reader.readLine()) != null) {
			stringBuilder.append(line);
			stringBuilder.append(ls);
		}
		
		reader.close();

		return stringBuilder.toString();
	}
	
	private String readFile(InputStream file) throws IOException {
		
		BufferedReader reader = new BufferedReader(new InputStreamReader(file));
		String line = null;
		StringBuilder stringBuilder = new StringBuilder();
		String ls = System.getProperty("line.separator");

		while ((line = reader.readLine()) != null) {
			stringBuilder.append(line);
			stringBuilder.append(ls);
		}
		
		reader.close();

		return stringBuilder.toString();
	}
}
