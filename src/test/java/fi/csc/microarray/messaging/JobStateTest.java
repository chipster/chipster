package fi.csc.microarray.messaging;

import static org.junit.Assert.*;

import org.junit.Test;

import fi.csc.microarray.messaging.JobState;

public class JobStateTest {

	@Test
	public void test() {
		
		assertTrue(JobState.COMPLETED.isFinished());
		assertTrue(JobState.FAILED.isFinished());
		assertTrue(JobState.FAILED_USER_ERROR.isFinished());
		assertTrue(JobState.ERROR.isFinished());
		assertTrue(JobState.CANCELLED.isFinished());
		assertTrue(JobState.TIMEOUT.isFinished());		
			
		assertFalse(JobState.NEW.isFinished());
		assertFalse(JobState.WAITING.isFinished());
		assertFalse(JobState.RUNNING.isFinished());
	}

}
