package fi.csc.chipster.web.adminweb.data;

import java.util.concurrent.locks.Lock;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import com.vaadin.ui.Label;
import com.vaadin.ui.Notification;
import com.vaadin.ui.Notification.Type;

import fi.csc.chipster.web.adminweb.ui.ReportView;

public class ReportDataSource {
	
	private static final Logger logger = Logger.getLogger(ReportDataSource.class);
	
	private StorageAdminAPI adminEndpoint;

	public ReportDataSource(StorageAdminAPI adminEndpoint) throws InstantiationException,
	IllegalAccessException {
		this.adminEndpoint = adminEndpoint;
	}
	
	public void update(final ReportView view) {
		
		String report;
		try {
			report = adminEndpoint.getStatusReport();		

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
			
		} catch (JMSException | InterruptedException e) {
			logger.error(e);
		}			
	}
}