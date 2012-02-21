package fi.csc.microarray.analyser.bsh;

import java.util.List;

import org.apache.log4j.Logger;

import bsh.EvalError;
import bsh.Interpreter;
import bsh.TargetError;
import fi.csc.microarray.analyser.ToolDescription;
import fi.csc.microarray.analyser.JobCancelledException;
import fi.csc.microarray.analyser.OnDiskAnalysisJobBase;
import fi.csc.microarray.analyser.ToolDescription.ParameterDescription;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.JobMessage.ParameterSecurityPolicy;
import fi.csc.microarray.messaging.message.JobMessage.ParameterValidityException;
import fi.csc.microarray.util.Exceptions;


/**
 * AnalysisJob for running BeanShell jobs.
 * 
 * A new BeanShell interpreter is instantiated for every job.
 * 
 * TODO Add better error handling. The bean shell script should be able to 
 * set job state, error message and output text
 * 
 * @author Taavi Hupponen, Aleksi Kallio
 * 
 */
public class BeanShellJob extends OnDiskAnalysisJobBase {

	public static class BSParameterSecurityPolicy implements ParameterSecurityPolicy {
		
		private static final int MAX_VALUE_LENGTH = 10000;
		
		public boolean isValueValid(String value, ParameterDescription parameterDescription) {
			
			// Check parameter size (DOS protection)
			if (value.length() > MAX_VALUE_LENGTH) {
				return false;
			}
			
			// No need to check content, parameters are passed inside Java Strings
			return true;
		}

	}
	
	public static BSParameterSecurityPolicy BS_PARAMETER_SECURITY_POLICY = new BSParameterSecurityPolicy();

	private static final Logger logger = Logger.getLogger(BeanShellJob.class);
	
	/**
	 * Wrap job information, create the BeanShell interpreter,
	 * pass the job info, and then run the script. 
	 * 
	 * 
	 */
	@Override
	protected void execute() throws JobCancelledException {
		updateStateDetailToClient("preparing BeanShell");
		
		// wrap the information to be passed to bean shell
		BeanShellJobInfo jobInfo = new BeanShellJobInfo();
		jobInfo.workDir = jobWorkDir;
		
		int i = 0;
		List<String> parameters;
		try {
			parameters = inputMessage.getParameters(BS_PARAMETER_SECURITY_POLICY, analysis);
		} catch (ParameterValidityException e) {
			outputMessage.setErrorMessage(e.getMessage()); // always has a message
			outputMessage.setOutputText(Exceptions.getStackTrace(e));
			updateState(JobState.FAILED_USER_ERROR, "");
			return;
		}
		for (ToolDescription.ParameterDescription param : analysis.getParameters()) {
			jobInfo.parameters.put(param.getName(), parameters.get(i));
			i++;
		}
		

		// create bean shell interpreter
		Interpreter interpreter = new Interpreter();
		
		try {
			// pass job info
			interpreter.set("jobInfo", jobInfo);
			
			// run the script
			updateStateDetailToClient("running the BeanShell script");
			interpreter.eval(analysis.getSourceCode());
		} 

		// analysis failed
		catch (TargetError te) {
			String errorMessage = "Running the BeanShell script failed.";
			logger.warn(errorMessage, te);
			outputMessage.setErrorMessage(errorMessage);
			outputMessage.setOutputText(Exceptions.getStackTrace(te));
			updateState(JobState.FAILED, "");
			return;
		} 
		
		// evaluation error
		catch (EvalError ee) {
			String errorMessage = "The BeanShell script could not be evaluated.";
			outputMessage.setErrorMessage(errorMessage);
			outputMessage.setOutputText(Exceptions.getStackTrace(ee));
			updateState(JobState.ERROR, "");
			return;
		}
		updateState(JobState.RUNNING, "BeanShell finished succesfully");
		
	}

	// TODO cancel by interrupting the interpreter somehow?
	@Override
	protected void cancelRequested() {
	}

}
