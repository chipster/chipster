/*
 * Created on Jan 27, 2005
 *
 */
package fi.csc.microarray.analyser.python;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.ToolDescription;
import fi.csc.microarray.analyser.JobCancelledException;
import fi.csc.microarray.analyser.OnDiskAnalysisJobBase;
import fi.csc.microarray.analyser.ProcessPool;
import fi.csc.microarray.analyser.ToolDescription.ParameterDescription;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.JobMessage.ParameterSecurityPolicy;
import fi.csc.microarray.messaging.message.JobMessage.ParameterValidityException;
import fi.csc.microarray.util.Exceptions;
import fi.csc.microarray.util.IOUtils;

/**
 * Uses Python to run actual analysis operations.
 * 
 * @author Aleksi Kallio
 */
public class PythonAnalysisJob extends OnDiskAnalysisJobBase {

	public static final String STRING_DELIMETER = "'";


	/**
	 * Checks that parameter values are safe to insert into Python code.
	 * Should closely match the code that is used to output the values in transformVariable(...).
	 * 
	 * @see PythonAnalysisJob#transformVariable(ParameterDescription, String)
	 *
	 */
	public static class PythonParameterSecurityPolicy implements ParameterSecurityPolicy {
		
		private static final int MAX_VALUE_LENGTH = 1000;
		
		/**
		 * This regular expression is very critical, because it checks code that is directly inserted
		 * into Python script. Hence it should be very conservative.
		 * 
		 * Interpretation: Maybe minus, zero or more digits, maybe point, zero or more digits.
		 */
		public static String NUMERIC_VALUE_PATTERN = "-?\\d*\\.?\\d*";
		
		/**
		 * This regular expression is not very critical, because text is inserted inside string constant in Python code.
		 * It should however always be combined with additional check that string terminator is not contained,
		 * because that way the string constant can be escaped. However values can be used in later
		 * points of the script in very different situations (filenames, etc.) and should be kept as simple as possible.
		 * 
		 * Interpretation: Only word characters and some special symbols are allowed.
		 */
		public static String TEXT_VALUE_PATTERN = "[\\w+\\-_:\\.,*() ]*";
		
		/**
		 *  @see ParameterSecurityPolicy#isValueValid(String, ParameterDescription)
		 */
		public boolean isValueValid(String value, ParameterDescription parameterDescription) {
			
			// Check parameter size (DOS protection)
			if (value.length() > MAX_VALUE_LENGTH) {
				return false;
			}
			
			// Check parameter content (Python injection protection)
			if (parameterDescription.isNumeric()) {
				
				// Numeric value must match the strictly specified pattern
				return value.matches(NUMERIC_VALUE_PATTERN);
				
			} else {
				
				// First check for string termination
				if (value.contains(STRING_DELIMETER)) {
					return false;
				}
				
				// Text value must still match specified pattern
				return value.matches(TEXT_VALUE_PATTERN);
			}
			
		}

	}
	
	public static PythonParameterSecurityPolicy PARAMETER_SECURITY_POLICY = new PythonParameterSecurityPolicy();
	
	static final Logger logger = Logger.getLogger(PythonAnalysisJob.class);
		
	private CountDownLatch waitPythonLatch = new CountDownLatch(1);
	
	// injected by handler at right after creation
	private ProcessPool processPool;
	private Process process;
	
	
	private class PythonProcessMonitor implements Runnable {

		public final String ERROR_MESSAGE_TOKEN = "Traceback";
		
		private ArrayList<String> outputLines;

		public void run() {
			
			logger.debug("Python process monitor started.");
			outputLines = new ArrayList<String>();
			
			BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
			
			boolean readMore = true;
			try {
				for (String line = reader.readLine(); readMore ; line = reader.readLine()) {
					
					// read end of stream --> error
					if (line == null || line.contains(SCRIPT_FAILED_STRING)) {
						updateState(JobState.FAILED, "Python script failed");
						readMore = false;
					} 
					
					// read script successful
					else if (line.contains(SCRIPT_SUCCESSFUL_STRING)) {
						updateState(JobState.COMPLETED, "Python script finished successfully");
						readMore = false;
					}
					
					// read normal output
					else {
						outputLines.add(line);
					}
				}
				
				// read the error message and chipster note
				if (getState() == JobState.FAILED) {

					// find the error token
					int errorLineNumber = -1;
					for (int i = outputLines.size(); i > 0 && errorLineNumber == -1; i--) {
						if (outputLines.get(i-1).startsWith(ERROR_MESSAGE_TOKEN)) {
							errorLineNumber = i-1;
						}
					}
				
					// get lines starting from the error token
					if (errorLineNumber != -1) {
						String errorMessage = "";
						errorMessage += outputLines.get(errorLineNumber);
						for (int i = errorLineNumber + 1; i < outputLines.size(); i++) {
							errorMessage += outputLines.get(i) + "\n";
						}
						errorMessage = errorMessage.substring(0, errorMessage.lastIndexOf("\n"));
						errorMessage = errorMessage.trim();
						
						// check for chipster note
						if (errorMessage.contains(CHIPSTER_NOTE_TOKEN)) {
							errorMessage = errorMessage.substring(errorMessage.indexOf(CHIPSTER_NOTE_TOKEN) + CHIPSTER_NOTE_TOKEN.length());
							errorMessage = errorMessage.trim();
							updateState(JobState.FAILED_USER_ERROR, "");
						}
						
						outputMessage.setErrorMessage(errorMessage);
					}
				}
				
			} catch (IOException e) {
				// also canceling the job leads here 
				logger.debug("error in monitoring Python process.");
				updateState(JobState.ERROR, "reading Python output failed.");
			}

			waitPythonLatch.countDown();
		}

		public String getOutput() {
			StringBuilder output = new StringBuilder();
			for (String line: outputLines) {
				output.append(line + "\n");
			}
			return output.toString();
		}
	
	}
	

	protected PythonAnalysisJob() {
	}

	
	
	/**
	 * Executes the analysis. 
	 * 
	 * @throws IOException
	 * @throws MicroarrayException 
	 * @throws InterruptedException 
	 */
	protected void execute() throws JobCancelledException {
		cancelCheck();
		updateStateDetailToClient("initialising Python");
		
		List<BufferedReader> inputReaders = new ArrayList<BufferedReader>();

		// load handler initialiser
		inputReaders.add(new BufferedReader(new StringReader(analysis.getInitialiser())));

		// load work dir initialiser
		logger.debug("job work dir: " + jobWorkDir.getPath());
		inputReaders.add(new BufferedReader(new StringReader("import os\nos.chdir('" + jobWorkDir.getName() + "')\n")));

		// load input parameters		
		int i = 0; 
		List<String> parameterValues;
		try {
			parameterValues = inputMessage.getParameters(PARAMETER_SECURITY_POLICY, analysis);
			
		} catch (ParameterValidityException e) {
			outputMessage.setErrorMessage(e.getMessage()); // always has a message
			outputMessage.setOutputText(Exceptions.getStackTrace(e));
			updateState(JobState.FAILED_USER_ERROR, "");
			return;
		}
		for (ToolDescription.ParameterDescription param : analysis.getParameters()) {
			String value = new String(parameterValues.get(i));
			String pythonSnippet = transformVariable(param, value);
			logger.debug("added parameter (" +  pythonSnippet + ")");
			inputReaders.add(new BufferedReader(new StringReader(pythonSnippet)));
			i++;
		}

		/*
		 * Enable importing of common scripts. Python won't know the dir of the
		 * script file, because it gets it from standard input.
		 */
		String setCommonsDir = 
				"import sys\n"
				+ "sys.path.append(chipster_common_path)\n";
		inputReaders.add(new BufferedReader(new StringReader(setCommonsDir)));
		
		// load input script
		String script = (String)analysis.getImplementation();
		inputReaders.add(new BufferedReader(new StringReader(script)));
		
		// load script finished trigger
		inputReaders.add(new BufferedReader(new StringReader("print '" + SCRIPT_SUCCESSFUL_STRING + "'\n")));
		
		
		// get a process
		cancelCheck();
		logger.debug("getting a process.");;
		try {
			this.process = processPool.getProcess();
		} catch (Exception e) {
			outputMessage.setErrorMessage("Starting Python failed.");
			outputMessage.setOutputText(Exceptions.getStackTrace(e));
			updateState(JobState.ERROR, "");
			return;
		}
		boolean processAlive = false;
		try {
			process.exitValue();
		} catch (IllegalThreadStateException itse) {
			processAlive = true;
		}
		if (!processAlive) {
			outputMessage.setErrorMessage("Starting Python failed.");
			String output = "";
			BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
			try {
				for (String line = reader.readLine(); line != null; line = reader.readLine()) {
					output += line + "\n";
				}
				reader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
				for (String line = reader.readLine(); line != null; line = reader.readLine()) {
					output += line + "\n";
				}
			} catch (IOException e) {
				logger.warn("could not read output stream");
			}
			outputMessage.setOutputText("Python already finished.\n\n" + output);
			updateState(JobState.ERROR, "");
			return;
		}
		
		updateStateDetailToClient("running Python");

		
		// launch the process monitor
		cancelCheck();
		logger.debug("about to start the Python process monitor.");
		PythonProcessMonitor processMonitor = new PythonProcessMonitor();
		new Thread(processMonitor).start();
		
		// combine the inputs into a single string and store it as the source code
		StringBuilder inputStringBuilder = new StringBuilder();
		try {
			for (BufferedReader reader : inputReaders) {			
				for (String line = reader.readLine(); line != null; line = reader.readLine()) {
					inputStringBuilder.append(line + "\n");
				}
			}
		} catch (IOException ioe) {
			logger.warn("creating Python input failed");
		}
	
		outputMessage.setSourceCode(inputStringBuilder.toString());
		
		// write the input to process
		logger.debug("writing the input to Python.");
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new OutputStreamWriter(process.getOutputStream()));
			writer.write(inputStringBuilder.toString());
			writer.newLine();
			writer.flush();
		} catch (IOException ioe) {
			// this happens if Python interpreter has died before or dies while writing the input
			// process monitor will notice this and set state etc
			logger.debug("writing input failed", ioe);
		} finally {
			IOUtils.closeIfPossible(writer);
		}
		
		// wait for the script to finish
		cancelCheck();
		logger.debug("waiting for the script to finish.");
		try {
			waitPythonLatch.await();
		} catch (InterruptedException e) {
			outputMessage.setErrorMessage("Running Python was interrupted.");
			outputMessage.setOutputText(Exceptions.getStackTrace(e));
			updateState(JobState.ERROR, "");
			return;
		}
		
		// script now finished or timeout
		cancelCheck();
		logger.debug("done waiting for " + analysis.getID() + ", state is " + getState());		

		// add output to result message
		String output = processMonitor.getOutput();
		outputMessage.setOutputText(output);
		
		// deal with errors and timeout if needed
		switch (getState()) {
		
		case RUNNING:
			outputMessage.setErrorMessage("Python did not finish before timeout.");
			updateState(JobState.TIMEOUT, "");
			return;
		case COMPLETED:
			// set state back to running, notify client
			updateState(JobState.RUNNING, "Python script finished successfully");
			updateStateDetailToClient("Python script finished successfully");
			return;
		case FAILED:
			// set error message if there is no specific message set already
			if (outputMessage.getErrorMessage() == null || outputMessage.getErrorMessage().equals("")) {
				outputMessage.setErrorMessage("Running Python script failed.");
			}
			return;
		case FAILED_USER_ERROR:
			return;
		case ERROR:
			// ProcessMonitor error
			outputMessage.setErrorMessage("Reading Python output failed.");
			return;
		default: 
			throw new IllegalStateException("Illegal job state: " + getState());
		}
	}

	@Override
	protected void preExecute() throws JobCancelledException {
		super.preExecute();
	}
	
	protected void postExecute() throws JobCancelledException {
		super.postExecute();
	}
	
	
	protected void cleanUp() {
		try {
			// release the process (don't recycle)
			if (process != null) {
				processPool.releaseProcess(process, false);
			}
		} catch (Exception e) {
			logger.error("error when releasing process. ", e);
		} finally {
			super.cleanUp();
		}
	}
	
	
	
	/**
	 * Converts a name-value -pair into Python variable definition.
	 */
	public static String transformVariable(ToolDescription.ParameterDescription param, String value) {
		
		// Escape strings and such
		if (!param.isNumeric()) {
			value = STRING_DELIMETER + value + STRING_DELIMETER; 
		}
		
		// If numeric, check for empty value
		if (param.isNumeric() && value.trim().isEmpty()) {
			value = "float('NaN')"; 
		}
		
		// Sanitize parameter name (remove spaces)
		String name = param.getName();
		name = name.replaceAll(" ", "_"); 
		
		// Construct and return parameter assignment
		return (name + " = " + value);
	}



	@Override
	protected void cancelRequested() {
		this.waitPythonLatch.countDown();
		
	}
	
	public void setProcessPool(ProcessPool processPool) {
		this.processPool = processPool;
	}
}
