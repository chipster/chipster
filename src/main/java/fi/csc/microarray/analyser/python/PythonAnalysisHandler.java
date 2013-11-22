package fi.csc.microarray.analyser.python;

import java.io.IOException;
import java.util.HashMap;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.InterpreterAnalysisHandler;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.analyser.ToolDescription;
import fi.csc.microarray.messaging.message.JobMessage;

/**
 * Handler for analysis tools written in Python.
 * 
 * @author Aleksi Kallio
 *
 */
public class PythonAnalysisHandler extends InterpreterAnalysisHandler {


	/**
	 * Logger for this class
	 */
	static final Logger logger = Logger
			.getLogger(PythonAnalysisHandler.class);

	

	public PythonAnalysisHandler(HashMap<String, String> parameters)
			throws IOException {
		super(parameters);
	}


	public AnalysisJob createAnalysisJob(JobMessage message, ToolDescription description, ResultCallback resultHandler) {
		PythonAnalysisJob analysisJob = new PythonAnalysisJob();
		analysisJob.construct(message, description, resultHandler);
		analysisJob.setProcessPool(this.processPool);
		return analysisJob;
	}

	@Override
	protected String getSingleCommentLineSymbol() {
		return "#";
	}


	@Override
	protected String getStringDelimeter() {
		return PythonAnalysisJob.STRING_DELIMETER;
	}


	@Override
	protected String getVariableNameSeparator() {
		return "_";
	}

}
