package fi.csc.microarray.analyser;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import org.testng.Assert;
import org.testng.annotations.AfterSuite;
import org.testng.annotations.BeforeSuite;

import fi.csc.microarray.TestConstants;
import fi.csc.microarray.client.tasks.TaskException;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.TaskEventListener;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.client.tasks.Task.State;
import fi.csc.microarray.databeans.fs.FSDataManager;
import fi.csc.microarray.messaging.MessagingTestBase;

/**
 * @author hupponen
 *
 */
public class AnalysisTestBase extends MessagingTestBase {

	protected FSDataManager manager;
	protected TaskExecutor executor;
	
	
	class JobResultListener implements TaskEventListener {

		private CountDownLatch latch;
		
		public JobResultListener(CountDownLatch latch) {
			this.latch = latch;
		}

		public void onStateChange(Task job, State oldState, State newState) {
			if (newState.isFinished()) {
				latch.countDown();
			}
		}
	}


	

	@BeforeSuite(alwaysRun = true)
	protected void setUp() throws Exception {
		super.setUp();
		this.manager = new FSDataManager();
		this.executor = new TaskExecutor(super.endpoint, manager);
		
	}

	@AfterSuite(alwaysRun = true)
	protected void tearDown() throws Exception {
		super.tearDown();
	}

	
	protected void executeJob(Task job) throws TaskException, InterruptedException {
		executeJob(job, TestConstants.TIMEOUT_AFTER, TimeUnit.MILLISECONDS, State.COMPLETED);
	}

	protected void executeJob(Task job, long timeout, TimeUnit timeUnit) throws TaskException, InterruptedException {
		executeJob(job, timeout, timeUnit, State.COMPLETED);
	}
	
	protected void executeJob(Task job, State expectedEndState) throws TaskException, InterruptedException {
		executeJob(job, TestConstants.TIMEOUT_AFTER, TimeUnit.MILLISECONDS, expectedEndState);
	}

	protected void executeJob(Task job, long timeout, TimeUnit timeUnit, State expectedEndState) throws TaskException, InterruptedException {
		CountDownLatch latch = new CountDownLatch(1);

		// send the job
		job.addTaskEventListener(new JobResultListener(latch));
		executor.startExecuting(job);
		
		// wait
		latch.await(timeout, timeUnit);

		// check end state
		State endState = job.getState();
		Assert.assertEquals(endState, expectedEndState);
	}

}
