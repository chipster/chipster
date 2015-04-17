package fi.csc.chipster.web.adminweb.data;

import java.io.IOException;
import java.io.Serializable;
import java.util.Collection;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import com.vaadin.data.util.BeanItemContainer;

import fi.csc.chipster.web.adminweb.ui.JobsView;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.admin.JobmanagerAdminAPI;
import fi.csc.microarray.messaging.admin.JobmanagerAdminAPI.JobsListener;
import fi.csc.microarray.messaging.admin.JobsEntry;

public class JobsContainer extends BeanItemContainer<JobsEntry> implements Serializable, JobsListener {
	
	private static final Logger logger = Logger.getLogger(JobsContainer.class);

	public static final String USERNAME = "username";
	public static final String OPERATION = "operation";
	public static final String STATUS = "status";
	public static final String COMPHOST = "compHost";
	public static final String START_TIME = "startTime";
	public static final String CANCEL_LINK = "cancelLink";

	public static final Object[] NATURAL_COL_ORDER  = new String[] {
		USERNAME, 		OPERATION, 		STATUS, 	COMPHOST, 		START_TIME, 	CANCEL_LINK };

	public static final String[] COL_HEADERS_ENGLISH = new String[] {
		"Username", 	"Operation", 	"Status", 	"Comp host", 	"Start time", 	"" };
	
	private JobmanagerAdminAPI jobmanagerAdminAPI;

	private JobsView view;




	public JobsContainer(JobsView view, JobmanagerAdminAPI jobmanagerAdminAPI) throws IOException, IllegalConfigurationException, MicroarrayException, JMSException {
		super(JobsEntry.class);
		this.view = view;
		this.jobmanagerAdminAPI = jobmanagerAdminAPI;
	}
	
	public void update() {		
		
		try {
			Collection<JobsEntry> list = jobmanagerAdminAPI.queryRunningJobs().values();
			statusUpdated(list);
						
		} catch (JMSException | InterruptedException | MicroarrayException e) {
			logger.error(e);
		}
	}

	@Override
	public void statusUpdated(final Collection<JobsEntry> jobs) {
		view.updateUI(new Runnable() {
			@Override
			public void run() {
				removeAllItems();
				
				for (JobsEntry entry : jobs) {
					addBean(entry);
				}
			}
		});
	}
}