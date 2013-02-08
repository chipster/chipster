package fi.csc.microarray.manager.web.ui;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.hibernate.Session;
import org.hibernate.exception.GenericJDBCException;

import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.server.ThemeResource;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.CheckBox;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.Panel;
import com.vaadin.ui.Table;
import com.vaadin.ui.VerticalLayout;

import fi.csc.microarray.manager.web.ChipsterAdminUI;
import fi.csc.microarray.manager.web.data.JobLogContainer;
import fi.csc.microarray.manager.web.data.JobLogEntry;
import fi.csc.microarray.manager.web.data.StatDataSource;
import fi.csc.microarray.manager.web.hbncontainer.HibernateUtil;


public class StatView extends VerticalLayout implements ClickListener {

	private Panel latestJobsPanel = new Panel("Latest jobs");
	private Panel jobCountsPanel = new Panel("Job counts monthly");
	private Panel topUsersPanel = new Panel("Top users during last year");
	private Panel toolFailsPanel = new Panel("Tool fails during last year");
	private Table latestJobsTable;
	private Table jobCountsTable;
	private Table topUsersTable;
	private Table toolFailsTable;

	private HorizontalLayout toolbarLayout;
	private HorizontalLayout statLayout = new HorizontalLayout();

	private ChipsterAdminUI app;

	private Button refreshButton;

	private CheckBox ignoreTestAccounts;

	private Session session;
	private StatDataSource dataSource;

	public StatView(ChipsterAdminUI app) {
		
		this.app = app;
					
		this.addComponent(getToolbar());		
		
		updateData(ignoreTestAccounts.getValue());
		
		latestJobsPanel.setHeight(100, Unit.PERCENTAGE);
		jobCountsPanel.setHeight(100, Unit.PERCENTAGE);
		topUsersPanel.setHeight(100, Unit.PERCENTAGE);
		toolFailsPanel.setHeight(100, Unit.PERCENTAGE);
		
		statLayout.addComponent(latestJobsPanel);
		statLayout.addComponent(jobCountsPanel);
		statLayout.addComponent(topUsersPanel);
		statLayout.addComponent(toolFailsPanel);
		
		this.addComponent(statLayout);
		
		this.setExpandRatio(statLayout, 1);
		
		this.setSizeFull();
	}
	
	private Session getHibernateSession() {
		if (session == null) {
			try {
				session = HibernateUtil.getSessionFactory().openSession();
			} catch (GenericJDBCException e) {
				//FIXME Show exception message and hide or disable all database based content
				e.printStackTrace();			
			}
		}
		return session;
	}

	private void updateData(boolean ignoreTestAccounts) {
		
		Session session = getHibernateSession();
		
		if (dataSource == null) {
			dataSource = new StatDataSource();
		}
		
		latestJobsTable = JobLogEntryToTable(dataSource.getLatestJobs(session));
		jobCountsTable = mapListToTable(dataSource.getJobsByMonth(session, ignoreTestAccounts));
		topUsersTable = mapListToTable(dataSource.getTopUsers(session, ignoreTestAccounts));
		toolFailsTable = mapListToTable(dataSource.getToolFails(session, ignoreTestAccounts));
		
		jobCountsTable.setVisibleColumns(dataSource.getJobsByMonthColumnOrder());
		topUsersTable.setVisibleColumns(dataSource.getTopUsersColumnOrder());
		toolFailsTable.setVisibleColumns(dataSource.getToolFailsColumnOrder());
		
		latestJobsPanel.setContent(latestJobsTable);
		jobCountsPanel.setContent(jobCountsTable);
		topUsersPanel.setContent(topUsersTable);
		toolFailsPanel.setContent(toolFailsTable);
	}
	
	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			
			toolbarLayout = new HorizontalLayout();
			
			refreshButton = new Button("Refresh");
			refreshButton.addClickListener((ClickListener)this);
			refreshButton.setIcon(new ThemeResource("../runo/icons/32/reload.png"));
			refreshButton.setEnabled(true);
			toolbarLayout.addComponent(refreshButton);
								
			ignoreTestAccounts = new CheckBox("Ignore test accounts", true);
			ignoreTestAccounts.addStyleName("toolbar-component");
			toolbarLayout.addComponent(ignoreTestAccounts);
						
			ignoreTestAccounts.addValueChangeListener(new ValueChangeListener() {

				@Override
				public void valueChange(ValueChangeEvent arg0) {
					updateData(ignoreTestAccounts.getValue());
				}
			});
			
			Label spaceEater = new Label(" ");
			toolbarLayout.addComponent(spaceEater);
			toolbarLayout.setExpandRatio(spaceEater, 1);

			toolbarLayout.addComponent(app.getTitle());	
			
			toolbarLayout.setWidth("100%");
			toolbarLayout.setStyleName("toolbar");
		}

		return toolbarLayout;
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
		if (event.getButton() == refreshButton) {			
			updateData(ignoreTestAccounts.getValue());
		}
	}
}
