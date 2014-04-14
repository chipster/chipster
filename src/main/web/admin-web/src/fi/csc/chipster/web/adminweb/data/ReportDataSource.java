package fi.csc.chipster.web.adminweb.data;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.locks.Lock;

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
import fi.csc.microarray.messaging.admin.CompAdminAPI;
import fi.csc.microarray.messaging.admin.ServerAdminAPI;
import fi.csc.microarray.messaging.admin.ServerAdminAPI.StatusReportListener;
import fi.csc.microarray.messaging.admin.StorageAdminAPI;
import fi.csc.microarray.messaging.message.ServerStatusMessage;

public class ReportDataSource {
	
	private static final Logger logger = Logger.getLogger(ReportDataSource.class);
	
	private StorageAdminAPI storageAdminAPI;
	private CompAdminAPI compAdminAPI;
	
	public void updateFilebrokerReport(final ReportView view) {
				
		String report;
		try {						
			
			report = getStorageAdminAPI().getStatusReport();		

			if (report != null) {
				Label label = view.getFilebrokerLabel();
				//Following is null if data loading in this thread
				//was faster than UI initialisation in another thread
				if (label.getUI() != null) {
					Lock labelLock = label.getUI().getSession().getLockInstance();
					labelLock.lock();
					try {
						label.setValue(report);

					} finally {
						labelLock.unlock();
					}
				}		
			} else {
				Notification.show("Timeout", "Chipster filebroker server doesn't respond", Type.ERROR_MESSAGE);
				logger.error("timeout while waiting status report");
			}
			
		} catch (JMSException | InterruptedException | IOException | IllegalConfigurationException | MicroarrayException e) {
			logger.error("failed to update storage status report", e);
		}			
	}
	
	private ServerAdminAPI getStorageAdminAPI() throws IOException, IllegalConfigurationException, MicroarrayException, JMSException {
		if (storageAdminAPI == null) {
			storageAdminAPI = new StorageAdminAPI();
		}
		return storageAdminAPI;
	}

	public void updateCompReport(final ReportView view) {
		
		try {	
			getCompAdminAPI().getStatusReports(new CompStatusReportListener(view, this));		
			
		} catch (JMSException | InterruptedException | IOException | IllegalConfigurationException | MicroarrayException e) {
			logger.error("failed to update comp status reports", e);
		}			
	}
	
	private CompAdminAPI getCompAdminAPI() throws IOException, IllegalConfigurationException, MicroarrayException, JMSException {
		if (compAdminAPI == null) {
			compAdminAPI = new CompAdminAPI();			
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

			VerticalLayout layout = view.getCompLayout();
			//Following is null if data loading in this thread
			//was faster than UI initialisation in another thread
			if (layout.getUI() != null) {
				Lock labelLock = layout.getUI().getSession().getLockInstance();
				labelLock.lock();
				try {
					layout.removeAllComponents();
					
					for (ServerStatusMessage serverStatus : statuses) {
					
						Label title = view.createReportLabel("Comp " + serverStatus.getHost() + "   ");
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

				} finally {
					labelLock.unlock();
				}
			}							
		}
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
			view.updateData();
		}
	}
}