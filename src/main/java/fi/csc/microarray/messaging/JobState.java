/*
 * Created on Feb 25, 2005
 *
 */
package fi.csc.microarray.messaging;

import java.util.Arrays;
import java.util.List;

/**
 * @author Taavi Hupponen
 *
 */
public enum JobState { 
	NEW, 
	RUNNING, 
	COMPLETED, 
	FAILED, // R failed, script error etc
	FAILED_USER_ERROR, // R failed, script error etc, we know that the user is to blame
	ERROR, // something went horribly wrong, unexpected exceptions etc
	TIMEOUT, 
	CANCELLED,
	COMP_BUSY,
	RESCHEDULED,
	SCHEDULED,
	WAITING,
	EXPIRED_WAITING; 

	static List<JobState> finished = Arrays.asList(
			COMPLETED, 
			FAILED, 
			FAILED_USER_ERROR,
			ERROR,
			CANCELLED,
			TIMEOUT,
			EXPIRED_WAITING);

	static List<JobState> finishedByComp = Arrays.asList(
			COMPLETED, 
			FAILED, 
			FAILED_USER_ERROR,
			ERROR,
			TIMEOUT);
	
	public boolean isFinished() {
		return finished.contains(this);
	}

	public boolean isFinishedByComp() {
		return finishedByComp.contains(this);
	}


}