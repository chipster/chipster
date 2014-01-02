package fi.csc.microarray.analyser.r;

import java.io.IOException;
import java.util.HashMap;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.InterpreterAnalysisHandler;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.analyser.ToolDescription;
import fi.csc.microarray.messaging.message.JobMessage;

public class RAnalysisHandler extends InterpreterAnalysisHandler {


	/**
	 * Logger for this class
	 */
	static final Logger logger = Logger
			.getLogger(RAnalysisHandler.class);


	public RAnalysisHandler(HashMap<String, String> parameters)
			throws IOException {
		super(parameters);
	}


	public AnalysisJob createAnalysisJob(JobMessage message, ToolDescription description, ResultCallback resultHandler) {
		RAnalysisJob analysisJob = new RAnalysisJob();
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
		return RAnalysisJob.STRING_DELIMETER;
	}


	@Override
	protected String getVariableNameSeparator() {
		return ".";
	}

}
