/*
 * Created on Jan 27, 2005
 *
 */
package fi.csc.microarray.analyser.r;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import org.apache.log4j.Logger;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.analyser.AnalysisDescription;
import fi.csc.microarray.analyser.AnalysisException;
import fi.csc.microarray.analyser.OnDiskAnalysisJobBase;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.messaging.JobState;

/**
 * Uses R to run actual analysis operations.
 * 
 * @author hupponen, akallio
 */
public class RAnalysisJob extends OnDiskAnalysisJobBase {

	static final Logger logger = Logger.getLogger(RAnalysisJob.class);
	
	private static String SCRIPT_SUCCESFUL_STRING = "nami-script-finished-succesfully";
	private static String SCRIPT_FAILED_STRING = "nami-script-finished-unsuccesfully";
	
	private int rTimeout;
	private CountDownLatch waitRLatch = new CountDownLatch(1);
	
	private Process process;
	
	
	private class RProcessMonitor implements Runnable {

		private StringBuffer output;

		public void run() {
			
			logger.debug("R process monitor started.");
			output = new StringBuffer("");
			
			BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
			
			boolean readMore = true;
			try {
				for (String line = reader.readLine(); readMore ; line = reader.readLine()) {
					
					// read end of stream --> error
					if (line == null || line.contains(SCRIPT_FAILED_STRING)) {
						logger.debug("R monitor read: " + line);
						updateState(JobState.FAILED, "R script failed", false);
						readMore = false;
					} 
					
					// read script successful
					// TODO better pattern matching 
					else if (line.contains(SCRIPT_SUCCESFUL_STRING)) {
						updateState(JobState.COMPLETED, "R script finished successfully", false);
						readMore = false;
					}
					
					// read normal output
					else {
						output.append(line + "\n");
					}
				}
			} catch (IOException e) {
				// also canceling the job leads here 
				logger.debug("Error in monitoring R process.");
				updateState(JobState.ERROR, "R monitor error", false);
			}

			waitRLatch.countDown();
		}

		public String getOutput() {
			return output.toString();
		}
	
	}
	

	protected RAnalysisJob() {
		this.rTimeout = DirectoryLayout.getInstance().getConfiguration().getInt("comp", "r-timeout");
	}

	
	
	/**
	 * Executes the analysis. 
	 * 
	 * @throws IOException
	 * @throws MicroarrayException 
	 * @throws InterruptedException 
	 */
	protected void execute() throws IOException, MicroarrayException, InterruptedException  {
		cancelCheck();
		updateStateDetail("initialising R", true);
		
		List<BufferedReader> inputReaders = new ArrayList<BufferedReader>();

		// load static initialiser
		inputReaders.add(new BufferedReader(new StringReader(AnalysisDescription.getStaticInitialiser())));
		
		// load work dir initialiser
		logger.debug("Job work dir: " + jobWorkDir.getPath());
		inputReaders.add(new BufferedReader(new StringReader("setwd(\"" + jobWorkDir.getName() + "\")\n")));
		
		
		// load input parameters		
		int i = 0; 
		for (AnalysisDescription.ParameterDescription param : analysis.getParameters()) {
			String value = new String(inputMessage.getParameters().get(i));
			String rSnippet = transformVariable(param.getName(), value, param.isNumeric());
			logger.debug("added parameter (" +  rSnippet + ")");
			inputReaders.add(new BufferedReader(new StringReader(rSnippet)));
			i++;
		}

		// load input script
		String script = (String)analysis.getImplementation();
		inputReaders.add(new BufferedReader(new StringReader(script)));
		
		// load script finished trigger
		inputReaders.add(new BufferedReader(new StringReader("print(\"" + SCRIPT_SUCCESFUL_STRING + "\")\n")));
		
		
		// get a process
		cancelCheck();
		logger.debug("Getting a process.");;
		this.process = resultHandler.getProcessPool().getProcess();
		updateStateDetail("running R", true);

		
		// launch the process monitor
		cancelCheck();
		logger.debug("About to start the R process monitor.");
		RProcessMonitor processMonitor = new RProcessMonitor();
		new Thread(processMonitor).start();
		
		// write the input to process
		logger.debug("Writing the input to R.");
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(process.getOutputStream()));
		for (BufferedReader reader : inputReaders) {			
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				writer.write(line);
				writer.newLine();
			}
		}
		writer.flush();

		
		// wait for the script to finish
		cancelCheck();
		logger.debug("Waiting for the script to finish.");
		waitRLatch.await(rTimeout, TimeUnit.SECONDS);
		
		
		// script now finished or timeout
		cancelCheck();
		logger.debug("Done waiting for " + analysis.getFullName() + ", state is " + getState());		

		// add output to result message
		String output = processMonitor.getOutput();
		// remove the first line (setwd(...))
		output = output.substring(output.indexOf(("\n")));
		outputMessage.setOutputText(output);
		
		// deal with errors and timeout if needed
		switch (getState()) {
		
		case RUNNING:
			// TODO add execution time this far to the message
			updateState(JobState.TIMEOUT, "R script exceeded timeout", false);
			throw new AnalysisException("Timeout occured before the analysis finished.");
		case FAILED:
			// state already updated
			
			// add STDERR to the output text
			StringWriter errorWriter = new StringWriter();
			BufferedReader errorReader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
			for (String line = errorReader.readLine(); line != null; line = errorReader.readLine()) {
				errorWriter.write(line + "\n");
			}
			outputMessage.setOutputText(output + "\n" + errorWriter.toString());

			throw new AnalysisException("Error occured when executing the analysis.");
		case COMPLETED:
			// set state back to running, notify client
			updateState(JobState.RUNNING, "R script finished successfully", true);
			break;			
		default: 
				String errorString = "Unexpected process end state: " + getState().toString();
				logger.error(this.inputMessage.getMessageID() + ": " + errorString);
				throw new AnalysisException(errorString);
		}
	}

	@Override
	protected void preExecute() throws Exception {
		super.preExecute();
	}
	
	protected void postExecute() throws Exception {
		super.postExecute();
	}
	
	
	protected void cleanUp() {
		try {
			// only try to recycle the process if the script finished succesfully
			// don't recycle at the moment
			if (process != null) {
				//this.resultHandler.getProcessPool().releaseProcess(process, getState().equals(JobState.SUCCESS));
				this.resultHandler.getProcessPool().releaseProcess(process, false);
			}
		} catch (Exception e) {
			logger.error("Error when releasing process. ", e);
		} finally {
			super.cleanUp();
		}

	}
	
	
	
	/**
	 * Converts a name-value -pair into R variable definition.
	 */
	public static String transformVariable(String name, String value, boolean isNumeric) {
		if (!isNumeric) {
			value = "\"" + value + "\""; // escape strings and such
		}
		name = name.replaceAll(" ", "_"); // remove spaces
		return (name + " <- " + value);
	}



	@Override
	protected void cancelRequested() {
		this.waitRLatch.countDown();
		
	}
}
