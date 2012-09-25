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

import fi.csc.microarray.manager.web.ChipsterAdminApplication;
import fi.csc.microarray.manager.web.data.JobLogContainer;
import fi.csc.microarray.manager.web.data.JobLogEntry;
import fi.csc.microarray.manager.web.data.StatDataSource;
import fi.csc.microarray.manager.web.hbncontainer.JobLogSessionManager;


public class StatView extends HorizontalLayout implements ClickListener {

	private ChipsterAdminApplication app;
	private Button testButton;

	private Panel latestJobsPanel = new Panel("Latest jobs");
	private Panel jobCountsPanel = new Panel("Job counts monthly");
	private Panel topUsersPanel = new Panel("Top users during last year");
	private Panel toolFailsPanel = new Panel("Tool fails during last year");
	private Table latestJobsTable;
	private Table jobCountsTable;
	private Table topUsersTable;
	private Table toolFailsTable;


	public StatView(ChipsterAdminApplication app) {
		this.app = app;

//		testButton = new Button("Test");
//		testButton.addListener(this);
//
//		this.addComponent(testButton);
		
		Session session = null;
		try {
			session = new JobLogSessionManager(app).getSession();
		} catch (GenericJDBCException e) {
			//FIXME Show exception message and hide or disable all database based content 
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
		
		latestJobsPanel.addComponent(latestJobsTable);
		jobCountsPanel.addComponent(jobCountsTable);
		topUsersPanel.addComponent(topUsersTable);
		toolFailsPanel.addComponent(toolFailsTable);
		
		latestJobsPanel.setHeight(100, UNITS_PERCENTAGE);
		jobCountsPanel.setHeight(100, UNITS_PERCENTAGE);
		topUsersPanel.setHeight(100, UNITS_PERCENTAGE);
		toolFailsPanel.setHeight(100, UNITS_PERCENTAGE);

		
		this.addComponent(latestJobsPanel);
		this.addComponent(jobCountsPanel);
		this.addComponent(topUsersPanel);
		this.addComponent(toolFailsPanel);
		
		this.setSizeFull();
	}

	private String mapListToHtml(List<Map> mapList) {

		StringBuffer html = new StringBuffer();

		if (mapList.size() > 0) {

			List<Object> keyList = new LinkedList<Object>(mapList.get(0).keySet());


			html.append("<table>");

			//Column headers
			html.append("<tr>");
			
			for (Object columnHeader : keyList) {
				html.append("<th>" + columnHeader + "</th>");
			}
			html.append("</tr>");

			
			//Content rows
			for (Map map : mapList) {
				html.append("<tr>");

				for (Object key : keyList) {
					html.append("<td>" + map.get(key) + "</td>");
				}
				html.append("</tr>");
			}


			html.append("</table>");
		}
		return html.toString();
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
			
				table.addItem(new Object[] { entry.getUsername(), entry.getOperation(), entry.getStartTime(), entry.getStatus() }, i++);
			}
		}
		return table;
	}
	
	private Table mapListToTable(List<Map> mapList) {

		Table table = new Table();
		table.setSizeFull();
		
		if (mapList.size() > 0) {

			List<Object> keyList = new LinkedList<Object>(mapList.get(0).keySet());


			//Column headers
			for (Object columnHeader : keyList) {
				table.addContainerProperty(columnHeader, String.class, null);
			}
			
			
			//Content rows
			int i = 0;
			for (Map map : mapList) {
			
				//We assume that the order of returned values is identical in all these maps 
				//otherwise values are put to wrong columns
				table.addItem(map.values().toArray(), i++);
			}
		}
		return table;
	}

	public void buttonClick(ClickEvent event) {
		if (event.getButton() == testButton) {

		}
	}
}
