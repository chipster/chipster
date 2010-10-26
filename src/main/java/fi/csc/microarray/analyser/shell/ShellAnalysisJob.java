package fi.csc.microarray.analyser.shell;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.util.LinkedList;
import java.util.concurrent.CountDownLatch;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.JobCancelledException;
import fi.csc.microarray.analyser.OnDiskAnalysisJobBase;
import fi.csc.microarray.analyser.AnalysisDescription.InputDescription;
import fi.csc.microarray.analyser.AnalysisDescription.OutputDescription;
import fi.csc.microarray.analyser.AnalysisDescription.ParameterDescription;
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
public class ShellAnalysisJob extends OnDiskAnalysisJobBase {
 
	protected String[] command;
	
    private Boolean useStdout = false;
    
    // Logger for this class
    static final Logger logger = Logger.getLogger(ShellAnalysisJob.class);
    
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
    	cancelCheck();
    	super.preExecute();
    	
        // Path to executable file
        String executablePath = analysis.getCommand();
        
        // Output method
        useStdout = analysis.getConfigParameters().get("stdout") != null &&
                analysis.getConfigParameters().get("stdout").toLowerCase().equals("yes");
        
        // Input parameter
        Boolean inputLast = analysis.getConfigParameters().get("input") != null &&
                analysis.getConfigParameters().get("input").toLowerCase().equals("last");
        
        // Additional arguments
        String arguments = analysis.getConfigParameters().get("arguments");
        String[] extraArguments = (arguments != null && !arguments.equals("")) ?
                arguments.split(",") : new String[] {};
 
        String outputParameter = null;
        if (!useStdout) {    
            // If program creates a normal file, we need an output parameter
            outputParameter = analysis.getConfigParameters().get("output");
        }

        LinkedList<String> inputParameters;
    	
    	// Get parameter values from user's input (order is significant)
        inputParameters = new LinkedList<String>(inputMessage.getParameters());
                
        // Generate the command to be executed
        LinkedList<String> commandParts = new LinkedList<String>();
        commandParts.add(executablePath);
        
        // Prepend arguments defined in the configuration file
        for (String arg : extraArguments) {
            commandParts.add(arg);
        }
        
        // Parameters
        int index = 0;
        for (ParameterDescription parameter : analysis.getParameters()) {
            String value = inputParameters.get(index);
            if (!value.equals("")) {
                commandParts.add("-" + parameter.getName());
                commandParts.add(value);
            }
            index++;
        }

        // Outputs to a file (currently we only support a single output)
        if (outputParameter != null) {
            OutputDescription output = analysis.getOutputFiles().get(0);
            commandParts.add("-" + outputParameter);
            commandParts.add(output.getFileName().getID());
        }
        
        // Inputs
        for (InputDescription input : analysis.getInputFiles()) {
        	if (!inputLast) {
                // Input is a named parameter
                commandParts.add("-" + input.getFileName());
            }
            commandParts.add(input.getFileName());
        }
        
        // store the command for execute()
        command = commandParts.toArray(new String[] {});
    }

    @Override
    protected void execute() throws JobCancelledException {
        try {
            String commandString = "";
        	for (String s : command) {
        		commandString += s;
        	}
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
            String outputString = Files.inputStreamToString(process.getInputStream()) + 
            					Files.inputStreamToString(process.getErrorStream());
            outputMessage.setOutputText(outputString);
            
            // failed job
            if (process.exitValue() != 0) {
                outputMessage.setErrorMessage("Running analysis tool failed.");
                updateState(JobState.FAILED, "non zero exit value");
                return;
            } 
        
            if (useStdout) {
                // FIXME not tested yet
                // Take output data from stdout
                InputStream stdoutStream =
                        new BufferedInputStream(process.getInputStream());
                OutputDescription output = analysis.getOutputFiles().get(0);
                File outputFile = new File(jobWorkDir, output.getFileName().getID());
                FileOutputStream fileStream = new FileOutputStream(outputFile);
                
                // Write from program's stdout to output file
                byte[] buffer = new byte[4096];  
                int bytesRead;  
                while ((bytesRead = stdoutStream.read(buffer)) != -1) {  
                    fileStream.write(buffer, 0, bytesRead);  
                }  
                stdoutStream.close();
                fileStream.close();
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
