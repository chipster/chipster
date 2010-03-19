package fi.csc.microarray.analyser.bsh;

import org.apache.log4j.Logger;

import bsh.EvalError;
import bsh.Interpreter;
import bsh.TargetError;
import fi.csc.microarray.analyser.AnalysisDescription;
import fi.csc.microarray.analyser.JobCancelledException;
import fi.csc.microarray.analyser.OnDiskAnalysisJobBase;
import fi.csc.microarray.messaging.JobState;


/**
 * AnalysisJob for running BeanShell jobs.
 * 
 * A new BeanShell interpreter is instantiated for every job.
 * 
 * @author hupponen
 * 
 */
public class BeanShellJob extends OnDiskAnalysisJobBase {

	private static final Logger logger = Logger.getLogger(BeanShellJob.class);
	
	/**
	 * Wrap job information, create the BeanShell interpreter,
	 * pass the job info, and then run the script. 
	 * 
	 * 
	 */
	@Override
	protected void execute() throws JobCancelledException {
		updateState(JobState.RUNNING, "preparing BeanShell", true);
		
		// wrap the information to be passed to bean shell
		BeanShellJobInfo jobInfo = new BeanShellJobInfo();
		jobInfo.workDir = jobWorkDir;
		
		int i = 0;
		for (AnalysisDescription.ParameterDescription param : analysis.getParameters()) {
			jobInfo.parameters.put(param.getName(), inputMessage.getParameters().get(i));
			i++;
		}
		

		// create bean shell interpreter
		Interpreter interpreter = new Interpreter();
		
		try {
			// pass job info
			interpreter.set("jobInfo", jobInfo);
			
			// run the script
			updateState(JobState.RUNNING, "running BeanShell", true);
			interpreter.eval(analysis.getSourceCode());
		} 

		// analysis failed
		catch (TargetError te) {
			String errorMessage = "Running the BeanShell script failed.";
			logger.warn(errorMessage, te);
			outputMessage.setErrorMessage(errorMessage);
			outputMessage.setOutputText(te.toString());
			updateState(JobState.FAILED, "", true);
			return;
		} 
		
		// evaluation error
		catch (EvalError ee) {
			String errorMessage = "The BeanShell script could not be evaluated.";
			outputMessage.setErrorMessage(errorMessage);
			outputMessage.setOutputText(ee.toString());
			updateState(JobState.ERROR, "", true);
			return;
		}
		updateState(JobState.RUNNING, "BeanShell finished succesfully", true);
		
	}

	// TODO cancel by interrupting the interpreter somehow?
	@Override
	protected void cancelRequested() {
	}

}
