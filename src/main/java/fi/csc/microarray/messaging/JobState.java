/*
 * Created on Feb 25, 2005
 *
 */
package fi.csc.microarray.messaging;


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
}
