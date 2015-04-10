package fi.csc.chipster.web.adminweb.data;

import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.Notification;
import com.vaadin.ui.Notification.Type;
import com.vaadin.ui.VerticalLayout;

import fi.csc.chipster.web.adminweb.ui.ReportView;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.admin.CompAdminAPI;
import fi.csc.microarray.messaging.admin.JobmanagerAdminAPI;
import fi.csc.microarray.messaging.admin.ServerAdminAPI;
import fi.csc.microarray.messaging.admin.ServerAdminAPI.StatusReportListener;
import fi.csc.microarray.messaging.admin.StorageAdminAPI;
import fi.csc.microarray.messaging.message.ServerStatusMessage;

public class ReportDataSource {
	
	private static final Logger logger = Logger.getLogger(ReportDataSource.class);
	
	private StorageAdminAPI storageAdminAPI;
	private CompAdminAPI compAdminAPI;
	private JobmanagerAdminAPI jobmanagerAdminAPI;

	private MessagingEndpoint endpoint;
	
	public ReportDataSource(MessagingEndpoint endpoint) {
		this.endpoint = endpoint;
	}

	public void updateFilebrokerReport(final ReportView view) {
				
		String report;
		try {						
			
			report = getStorageAdminAPI().getStatusReport();		

			if (report != null) {
				updateFilebrokerUI(view, report);
			} else {
				Notification.show("Timeout", "Chipster filebroker server doesn't respond", Type.ERROR_MESSAGE);
				logger.error("timeout while waiting status report");
			}
			
		} catch (JMSException | InterruptedException | IOException | IllegalConfigurationException | MicroarrayException e) {
			logger.error("failed to update storage status report", e);
		}			
	}
	
	private void updateFilebrokerUI(final ReportView view, final String report) {
		view.updateUI(new Runnable() {
			@Override
			public void run() {
				Label label = view.getFilebrokerLabel();
				label.setValue(report);						
			}
		});
	}

	private ServerAdminAPI getStorageAdminAPI() throws IOException, IllegalConfigurationException, MicroarrayException, JMSException {
		if (storageAdminAPI == null) {
			storageAdminAPI = new StorageAdminAPI(endpoint);
		}
		return storageAdminAPI;
	}
	
	private JobmanagerAdminAPI getJobmanagerAdminAPI() throws IOException, IllegalConfigurationException, MicroarrayException, JMSException {
		if (jobmanagerAdminAPI == null) {
			jobmanagerAdminAPI = new JobmanagerAdminAPI(endpoint);
		}
		return jobmanagerAdminAPI;
	}

	public void updateCompReport(final ReportView view, int timeout) {
		
		try {	
			getCompAdminAPI().getStatusReports(new CompStatusReportListener(view, this), timeout);		
			
		} catch (JMSException | InterruptedException | IOException | IllegalConfigurationException | MicroarrayException e) {
			logger.error("failed to update comp status reports", e);
		}			
	}
	
	public void updateJobmanagerReport(final ReportView view) {
		
		String report;
		try {						
			
			report = getJobmanagerAdminAPI().getStatusReport();		

			if (report != null) {
				
				updateJobmanagerUI(view, report);
				
			} else {
				Notification.show("Timeout", "Chipster jobmanager server doesn't respond", Type.ERROR_MESSAGE);
				logger.error("timeout while waiting status report");
			}
			
		} catch (JMSException | IOException | IllegalConfigurationException | MicroarrayException | InterruptedException e) {
			logger.error("failed to update storage status report", e);
		}		
	}
	
	private void updateJobmanagerUI(final ReportView view, final String report) {
		view.updateUI(new Runnable() {
			@Override
			public void run() {
				VerticalLayout layout = view.getJobmanagerLayout();

				layout.removeAllComponents();

				Button purgeButton = view.createReportButton("Purge old jobs");

				purgeButton.addClickListener(new PurgeClickListener(view, ReportDataSource.this));

				Label reportLabel = view.createReportLabel(report);
				layout.addComponent(reportLabel);
				layout.addComponent(purgeButton);				
			}
		});
	}

	private CompAdminAPI getCompAdminAPI() throws IOException, IllegalConfigurationException, MicroarrayException, JMSException {
		if (compAdminAPI == null) {
			compAdminAPI = new CompAdminAPI(endpoint);			
		}
		return compAdminAPI;
	}

	public static class CompStatusReportListener implements StatusReportListener {
		private ReportView view;
		private ReportDataSource reportDataSource;

		public CompStatusReportListener(ReportView view, ReportDataSource reportDataSource) {
			this.view = view;
			this.reportDataSource = reportDataSource;
		}

		@Override
		public void statusUpdated(List<ServerStatusMessage> statuses) {

			updateCompUI(view, statuses, reportDataSource);
		}
	}
	
	private static void updateCompUI(final ReportView view,
			final List<ServerStatusMessage> statuses, final ReportDataSource reportDataSource) {

		view.updateUI(new Runnable() {
			@Override
			public void run() {
				VerticalLayout layout = view.getCompLayout();

				layout.removeAllComponents();

				Collections.sort(statuses, new Comparator<ServerStatusMessage>() {
					@Override
					public int compare(ServerStatusMessage m1, ServerStatusMessage m2) {
						int hostComparison = m1.getHost().compareTo(m2.getHost());
						int idComparison = m1.getHostId().compareTo(m2.getHostId());
						if (hostComparison != 0) {
							return hostComparison;
						}
						return idComparison;
					}						
				});

				for (ServerStatusMessage serverStatus : statuses) {

					Label title = view.createReportLabel("Comp " + serverStatus.getHost() + " (" + serverStatus.getHostId() + ")" + "    ");
					Button shutdownButton = view.createReportButton("Stop gracefully");

					shutdownButton.addClickListener(new StopClickListener(view, reportDataSource, serverStatus.getHostId()));
					HorizontalLayout titleRow = new HorizontalLayout();						
					titleRow.addComponent(title);
					titleRow.addComponent(shutdownButton);
					Label reportLabel = view.createReportLabel(
							serverStatus.toString());
					layout.addComponent(titleRow);
					layout.addComponent(reportLabel);
				}
			}
		});
	}

	public static class StopClickListener implements ClickListener {

		private ReportView view;
		private ReportDataSource reportDataSource;
		private String hostId;

		public StopClickListener(ReportView view, ReportDataSource reportDataSource, String hostId) {
			this.view = view;
			this.reportDataSource = reportDataSource;
			this.hostId = hostId;
		}

		@Override
		public void buttonClick(ClickEvent event) {
			try {
				reportDataSource.getCompAdminAPI().stopGracefullyComp(hostId);
			} catch (IOException | IllegalConfigurationException | MicroarrayException | JMSException e) {
				e.printStackTrace();
			}
			view.update();
		}
	}
	
	public static class PurgeClickListener implements ClickListener {

		private ReportView view;
		private ReportDataSource reportDataSource;

		public PurgeClickListener(ReportView view, ReportDataSource reportDataSource) {
			this.view = view;
			this.reportDataSource = reportDataSource;
		}

		@Override
		public void buttonClick(ClickEvent event) {
			try {
				reportDataSource.getJobmanagerAdminAPI().purge();
			} catch (IOException | IllegalConfigurationException | MicroarrayException | JMSException e) {
				e.printStackTrace();
			}
			view.update();
		}
	}
}