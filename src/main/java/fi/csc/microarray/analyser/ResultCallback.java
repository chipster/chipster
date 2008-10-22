/*
 * Created on Mar 4, 2005
 *
 */

package fi.csc.microarray.analyser;

import java.io.File;

import fi.csc.microarray.messaging.message.NamiMessage;
import fi.csc.microarray.messaging.message.ResultMessage;

/**
 * For reporting back analysis results. These methods are called by the thread running a job.
 * 
 * @author hupponen
 */
public interface ResultCallback {

	
	/**
	 * All the data must be sent before returning, as after calling this method, the thread
	 * running the job will clean up job data files.
	 * 
	 * 
	 * @param inputMessage
	 * @param resultMessage
	 */
	public void sendResultMessage(NamiMessage inputMessage, ResultMessage resultMessage); 
	
	public File getWorkDir();
	public boolean shouldSweepWorkDir();
	
	public ProcessPool getProcessPool();
	
	public void removeRunningJob(AnalysisJob job);
}