package fi.csc.microarray.messaging;


/**
 *
 * @author  Aleksi Kallio
 */

public class AuthenticationTest extends MessagingTestBase {
	
//	private boolean gotData = false;
//	private TaskExecutor executor;
//	private DataBean microarray;
//	private DataManager manager; 
//	
//	public AuthenticationTest(String s) throws IOException {
//		this.manager = new DataManager();
//	}    
//
//	@BeforeTest
//	protected void setUp() throws Exception {
//		super.setUp();
//		microarray = manager.createDataBean("microarray1.tsv", this.getClass().getResourceAsStream("/microarray.tsv"));
//		this.executor = new TaskExecutor(super.endpoint, manager);
//	}
//	
//	@AfterTest
//	protected void tearDown() throws Exception {
//		super.tearDown();
//	}
//	
//	@Test
//	public void testAuthenticationAndSessions() throws JMSException, IOException, MicroarrayException, InterruptedException, TaskException{
//
//		// start job
//		Task job = executor.createTask(TestConstants.DUMMY_ANALYSIS_NAME);
//		job.addInput("microarray.tsv", microarray);
//		job.addTaskEventListener(new TaskEventListener() {
//			public void onStateChange(Task task, State oldState, State newState) {
//				if (newState == State.COMPLETED) {
//					gotData = true;
//				}
//			}
//
//		});
//		executor.startExecuting(job, TestConstants.TIMEOUT_AFTER);
//		
//		// wait for auth request
//		for (int i = 0; i < 4; i++) {
//			Thread.sleep(1024);
//			if (authenticationListener.isRequested()) {
//				continue;
//			}
//		}
//		
//		Assert.assertTrue(authenticationListener.isRequested());
//
//		// wait for data
//		for (int i = 0; i < 10; i++) {
//			Thread.sleep(1024);
//			if (gotData) {
//				continue;
//			}
//		}
//		Assert.assertTrue(gotData);
//		
//		// send another request
//		authenticationListener.reset();
//		gotData = false;
//		Task job2 = executor.createTask(TestConstants.DUMMY_ANALYSIS_NAME);
//		job2.addInput("microarray.tsv", microarray);
//		job2.addTaskEventListener(new TaskEventListener() {
//			public void onStateChange(Task task, State oldState, State newState) {
//				if (newState == State.COMPLETED) {
//					gotData = true;
//				}
//			}
//		});
//		executor.startExecuting(job2, TestConstants.TIMEOUT_AFTER);
//
//		// wait for auth request (which should not come this time)
//		for (int i = 0; i < 4; i++) {
//			Thread.sleep(1024);
//			if (authenticationListener.isRequested()) {
//				Assert.fail("authentication was requested but session should have existed");
//			}
//		}
//		Assert.assertTrue(!authenticationListener.isRequested());
//
//		// wait for data, again
//		for (int i = 0; i < 10; i++) {
//			Thread.sleep(1024);
//			if (gotData) {
//				continue;
//			}
//		}
//		Assert.assertTrue(gotData);
//	}	
}