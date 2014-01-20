package fi.csc.chipster.web.adminweb.data;

import java.io.Serializable;
import java.util.Date;

import com.vaadin.data.util.BeanItemContainer;

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
	
//	private static final String[] status = new String[] { "RUNNING", "WAITING", "TRANSFERRING FILES" }; 

	public void update() {
		
      JobsEntry entry = new JobsEntry();
      
      entry.setUsername("Not implemented yet");
      entry.setOperation("");
      entry.setStatus("");
      entry.setCompHost("");
      entry.setStartTime(new Date());
      
      this.addBean(entry);
		
//		Random rnd = new Random();
//		
//		int count = rnd.nextInt(10);
//		
//		
//		removeAllItems();
//		
//    	
//		JobsEntry entry;
//
//        for (int i = 0; i < count; i++) {
//            entry = new JobsEntry();
//            
//            entry.setUsername(RandomUtil.getRandomUserName(rnd));
//            entry.setOperation(RandomUtil.getRandomOperation(rnd));
//            entry.setStatus(status[rnd.nextInt(status.length)]);
//            entry.setCompHost(RandomUtil.getRandomComp(rnd));
//            entry.setStartTime(RandomUtil.getRandomDateToday(rnd));
//            
//            this.addBean(entry);
//        }
    }
}