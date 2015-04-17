package fi.csc.chipster.web.adminweb.ui;

import java.io.IOException;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;

import fi.csc.chipster.web.adminweb.ChipsterAdminUI;
import fi.csc.chipster.web.adminweb.data.JobsContainer;
import fi.csc.chipster.web.adminweb.util.Notificationutil;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.admin.JobmanagerAdminAPI;
import fi.csc.microarray.messaging.admin.JobsEntry;

public class JobsView extends AsynchronousView implements ClickListener {
	
	private static final Logger logger = Logger.getLogger(JobsView.class);
	
	public static final int WAIT_SECONDS = 5;
	
	private HorizontalLayout toolbarLayout;

	private JobsTable table;
	private JobmanagerAdminAPI jobmanagerAdminAPI;
	private JobsContainer dataSource;

	public JobsView(ChipsterAdminUI app) {
		
		super(app, WAIT_SECONDS);
		
		try {
			jobmanagerAdminAPI = new JobmanagerAdminAPI(app.getEndpoint());
			dataSource = new JobsContainer(this, jobmanagerAdminAPI);

			table = new JobsTable(this);
			table.setContainerDataSource(dataSource);

			table.setVisibleColumns(JobsContainer.NATURAL_COL_ORDER);
			table.setColumnHeaders(JobsContainer.COL_HEADERS_ENGLISH);

			this.addComponent(getToolbar());
			this.addComponent(super.getProggressIndicator());
			this.addComponent(table);

			setSizeFull();
			this.setExpandRatio(table, 1);

		} catch (IOException | IllegalConfigurationException | MicroarrayException | JMSException e) {
			logger.error("can't initialize jobs view", e);
		} 
	}

	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			
			toolbarLayout = new HorizontalLayout();
			
			toolbarLayout.addComponent(super.createRefreshButton(this));
			
			Label spaceEater = new Label(" ");
			toolbarLayout.addComponent(spaceEater);
			toolbarLayout.setExpandRatio(spaceEater, 1);
			
			toolbarLayout.addComponent(super.getApp().getTitle());	
			
			toolbarLayout.setWidth("100%");
			toolbarLayout.setStyleName("toolbar");
		}

		return toolbarLayout;
	}

	public void buttonClick(ClickEvent event) {
		final Button source = event.getButton();

		if (super.isRefreshButton(source)) {
			update();
		} 
	}
	
	public void update() {
		
		super.submitUpdate(new Runnable() {

			@Override
			public void run() {				
				dataSource.update();				
			}			
		});				
	}

	public void cancel(JobsEntry job) {
		try {
			jobmanagerAdminAPI.cancelJob(job.getJobId());
			table.removeItem(job);
		} catch (MicroarrayException e) {
			Notificationutil.showFailNotification(e.getClass().getSimpleName(), e.getMessage());
		}
	}

	public JobsTable getEntryTable() {
		return table;
	}
}
