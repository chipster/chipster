package fi.csc.microarray.manager.web.data;

import java.io.Serializable;
import java.util.Random;

import com.vaadin.data.util.BeanItemContainer;

import fi.csc.microarray.manager.web.util.RandomUtil;

public class JobsContainer extends BeanItemContainer<JobsEntry> implements
Serializable {

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



	public JobsContainer() {
		super(JobsEntry.class);
	}
	
	private static final String[] status = new String[] { "RUNNING", "WAITING", "TRANSFERING FILES" }; 

	public void update() {
		
		Random rnd = new Random();
		
		int count = rnd.nextInt(10);
		
		removeAllItems();
		
    	
		JobsEntry entry;

        for (int i = 0; i < count; i++) {
            entry = new JobsEntry();
            
            entry.setUsername(RandomUtil.getRandomUserName(rnd));
            entry.setOperation(RandomUtil.getRandomOperation(rnd));
            entry.setStatus(status[rnd.nextInt(status.length)]);
            entry.setCompHost(RandomUtil.getRandomComp(rnd));
            entry.setStartTime(RandomUtil.getRandomDateToday(rnd));
            
            this.addBean(entry);
        }
    }
}