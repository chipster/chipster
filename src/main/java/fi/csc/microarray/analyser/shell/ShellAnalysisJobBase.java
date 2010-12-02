package fi.csc.microarray.analyser.shell;

import java.io.File;
import java.io.FileOutputStream;
import java.util.concurrent.CountDownLatch;

import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.JobCancelledException;
import fi.csc.microarray.analyser.OnDiskAnalysisJobBase;
import fi.csc.microarray.analyser.AnalysisDescription.OutputDescription;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.util.Files;

/**
 * Job that is run as a generic shell command.
 * <p>
 * Parameters that could be given in the tool configuration:
 *  <ul>
 *  <li> output - name of the parameter for passing output files
 *  <li> stdout - if equals to "yes", the output is read from stdout
 *  <li> input - if equals to "last", then input is given without parameter
 *               name and as the last parameter
 *  <li> arguments - comma-separated list of arguments that will get passed
 *               as the first arguments for the command.   
 *  </ul>
 * 
 * @author naktinis, hupponen
 *
 */
public abstract class ShellAnalysisJobBase extends OnDiskAnalysisJobBase {
 
	protected String[] command;
	
    protected Boolean useStdout = false;
    
    // Logger for this class
    static final Logger logger = Logger.getLogger(ShellAnalysisJobBase.class);
    
    // Latch for canceling or finishing a job
    private CountDownLatch latch = new CountDownLatch(1);
    
    // Operating system process
    private Process process = null;

    
    /**
     * Construct the command line.
     * 
     */
    @Override
    protected void preExecute() throws JobCancelledException {
    	super.preExecute();
    }

    @Override
    protected void execute() throws JobCancelledException {
        try {
            String commandString = "";
        	for (String s : command) {
        		commandString += s + " ";
        	}
        	commandString = commandString.trim();
        	logger.info("running shell job: " + commandString);
            
            process = Runtime.getRuntime().exec(command, null, jobWorkDir);
            updateStateDetailToClient("running analysis tool");
            
            // Start a new thread to listen to OS process status
            new Thread(new ProcessWaiter()).start();
            
            // wait for the job to finish
            latch.await();
            cancelCheck();
            
            // now finished
            updateStateDetailToClient("analysis tool finished");
            
            // put stdout and stderr to outputmessage
            // stdoOutString is also possible needed later for the out file
            String stdOutString = Files.inputStreamToString(process.getInputStream());
            String stdErrorString = Files.inputStreamToString(process.getErrorStream());
            outputMessage.setOutputText(stdOutString + "\n" + stdErrorString);
            
            // failed job
            if (process.exitValue() != 0) {
                outputMessage.setErrorMessage("Running analysis tool failed.");
                updateState(JobState.FAILED, "non zero exit value");
                return;
            } 
        
            if (useStdout) {
            	// use screen output as result data
                OutputDescription output = analysis.getOutputFiles().get(0);
                File outputFile = new File(jobWorkDir, output.getFileName().getID());
                FileOutputStream fileStream = null;
                try {
                	fileStream = new FileOutputStream(outputFile);
                	IOUtils.write(stdOutString, fileStream);
                	fileStream.flush();
                } finally {
                	IOUtils.closeQuietly(fileStream);
                }
            }

            // if successful, don't need to do anything, just leave the state as running

        } catch (Exception e) {
        	outputMessage.setErrorMessage("Running analysis tool failed.");
        	outputMessage.setOutputText(e.toString());
        	updateState(JobState.ERROR, "analysis tool failed");
        	return;
        }
    }

    
    
    
    /**
     * User decided to cancel this job.
     */
    @Override
    protected void cancelRequested() {
        latch.countDown();
    }
    
    /**
     * Destroy operating system process if it is still
     * running.
     */
    @Override
    protected void cleanUp() {
		
    	// kill the process if not already finished
    	try {
			if (process != null) {
				try {
					process.exitValue();
				} catch (IllegalThreadStateException itse) {
					process.destroy();
				}
			}
		} catch (Exception e) {
			logger.error("error when destroying process ", e);
		} finally {
			super.cleanUp();
		}
    }


    /**
     * A simple runnable that waits for an operating system
     * process to finish and reduces a given latch by one.
     * 
     */
    private class ProcessWaiter implements Runnable {

		@Override
		public void run() {
    		try {
                process.waitFor();
            } catch (InterruptedException e) {
            	throw new RuntimeException(e);
            } finally {
            	latch.countDown();
            }
		}
    }
}
