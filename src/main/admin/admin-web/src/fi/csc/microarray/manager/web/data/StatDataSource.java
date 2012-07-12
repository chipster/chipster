package fi.csc.microarray.manager.web.data;

import java.util.Calendar;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.List;
import java.util.Map;

import org.hibernate.Session;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.transform.Transformers;

public class StatDataSource {
	
	public static final String ROW_COUNT = "rowCount";

	
	public Collection<JobLogEntryRowCount> getTopUsers(Session session) {
		
		//Start one year ago
		Calendar fromDate = new GregorianCalendar();
		fromDate.set(Calendar.YEAR, fromDate.get(Calendar.YEAR) - 1);
		
		@SuppressWarnings("unchecked")
		List<JobLogEntryRowCount> results = session.createCriteria(JobLogEntry.class)
				.add(Restrictions.ge(JobLogContainer.START_TIME, fromDate.getTime()))
			    .setProjection( Projections.projectionList()
			        .add(Projections.rowCount(), ROW_COUNT)
			        .add(Projections.groupProperty(JobLogContainer.USERNAME), JobLogContainer.USERNAME)		        
			    )
			    .addOrder(Order.desc(ROW_COUNT))
			    .setMaxResults(10)
			    .setResultTransformer(Transformers.aliasToBean(JobLogEntryRowCount.class))
			    .list();
		
		return results;
	}
	
	public Collection<JobLogEntry> getLastJobs(Session session) {
		
		@SuppressWarnings("unchecked")
		List<JobLogEntry> results = session.createCriteria(JobLogEntry.class)
			    .addOrder(Order.desc(JobLogContainer.START_TIME))
			    .setMaxResults(10)
			    .setResultTransformer(Transformers.aliasToBean(JobLogEntry.class))
			    .list();
		
		return results;
	}
	
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public List<Map> getJobsByMonth(Session session) {
		
		return session.createQuery("select new map(year(startTime) as year, month(startTime) as month, count(*) as count) " +
				"from JobLogEntry " + 
						"group by year(startTime), month(startTime) " +
						"order by year(startTime) asc, month(startTime) asc").list();
	}
}
