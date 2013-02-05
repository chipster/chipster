package fi.csc.microarray.manager.web.ui;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.hibernate.Session;
import org.hibernate.exception.GenericJDBCException;

import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Panel;
import com.vaadin.ui.Table;

import fi.csc.microarray.manager.web.data.JobLogContainer;
import fi.csc.microarray.manager.web.data.JobLogEntry;
import fi.csc.microarray.manager.web.data.StatDataSource;
import fi.csc.microarray.manager.web.hbncontainer.HibernateUtil;


public class StatView extends HorizontalLayout implements ClickListener {

	private Button testButton;

	private Panel latestJobsPanel = new Panel("Latest jobs");
	private Panel jobCountsPanel = new Panel("Job counts monthly");
	private Panel topUsersPanel = new Panel("Top users during last year");
	private Panel toolFailsPanel = new Panel("Tool fails during last year");
	private Table latestJobsTable;
	private Table jobCountsTable;
	private Table topUsersTable;
	private Table toolFailsTable;


	public StatView() {
		
		Session session = null;
		try {
			session = HibernateUtil.getSessionFactory().openSession();
		} catch (GenericJDBCException e) {
			//FIXME Show exception message and hide or disable all database based content
			e.printStackTrace();
			return;
		}

		StatDataSource dataSource = new StatDataSource();

		
		latestJobsTable = JobLogEntryToTable(dataSource.getLatestJobs(session));
		jobCountsTable = mapListToTable(dataSource.getJobsByMonth(session));
		topUsersTable = mapListToTable(dataSource.getTopUsers(session));
		toolFailsTable = mapListToTable(dataSource.getToolFails(session));
		
		jobCountsTable.setVisibleColumns(dataSource.getJobsByMonthColumnOrder());
		topUsersTable.setVisibleColumns(dataSource.getTopUsersColumnOrder());
		toolFailsTable.setVisibleColumns(dataSource.getToolFailsColumnOrder());
		
		latestJobsPanel.setContent(latestJobsTable);
		jobCountsPanel.setContent(jobCountsTable);
		topUsersPanel.setContent(topUsersTable);
		toolFailsPanel.setContent(toolFailsTable);
		
		latestJobsPanel.setHeight(100, Unit.PERCENTAGE);
		jobCountsPanel.setHeight(100, Unit.PERCENTAGE);
		topUsersPanel.setHeight(100, Unit.PERCENTAGE);
		toolFailsPanel.setHeight(100, Unit.PERCENTAGE);

		
		this.addComponent(latestJobsPanel);
		this.addComponent(jobCountsPanel);
		this.addComponent(topUsersPanel);
		this.addComponent(toolFailsPanel);
		
		this.setSizeFull();
	}
	
	private Table JobLogEntryToTable(List<JobLogEntry> list) {

		Table table = new Table();
		table.setSizeFull();
		
		if (list.size() > 0) {

			List<Object> keyList = Arrays.asList(new Object[] { JobLogContainer.USERNAME, JobLogContainer.OPERATION, JobLogContainer.START_TIME, JobLogContainer.STATUS });


			//Column headers
			for (Object columnHeader : keyList) {
				table.addContainerProperty(columnHeader, String.class, null);
			}
			
			
			//Content rows
			int i = 0;
			for (JobLogEntry entry : list) {
			
				table.addItem(new Object[] { entry.getUsername(), entry.getOperation(), entry.getStartTime().toString(), entry.getStatus() }, i++);
			}
		}
		return table;
	}
	
	private Table mapListToTable(List<Map<Object, Object>> list) {

		Table table = new Table();
		table.setSizeFull();
		
		if (list.size() > 0) {

			List<Object> keyList = new LinkedList<Object>(list.get(0).keySet());


			//Column headers
			for (Object columnHeader : keyList) {
				table.addContainerProperty(columnHeader, String.class, null);
			}
			
			
			//Content rows
			int i = 0;
			for (Map<Object, Object> map : list) {
			
				//We assume that the order of returned values is identical in all these maps 
				//otherwise values are put to wrong columns
				List<String> stringValues = new LinkedList<String>();
				
				for (Object objValue : map.values()) {
					if (objValue != null) {
						stringValues.add(objValue.toString());
					} else {
						stringValues.add("");
					}					
				}
				
				table.addItem(stringValues.toArray(), i++);
			}
		}
		return table;
	}

	public void buttonClick(ClickEvent event) {
		if (event.getButton() == testButton) {

		}
	}
}
