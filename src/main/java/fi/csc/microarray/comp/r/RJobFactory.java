package fi.csc.microarray.comp.r;

import java.io.IOException;
import java.util.HashMap;

import org.apache.log4j.Logger;

import fi.csc.chipster.toolbox.ToolboxTool;
import fi.csc.microarray.comp.CompException;
import fi.csc.microarray.comp.CompJob;
import fi.csc.microarray.comp.InterpreterJobFactory;
import fi.csc.microarray.comp.ResultCallback;
import fi.csc.microarray.comp.ToolDescription;
import fi.csc.microarray.messaging.message.GenericJobMessage;

public class RJobFactory extends InterpreterJobFactory {

	static final Logger logger = Logger
			.getLogger(RJobFactory.class);


	public RJobFactory(HashMap<String, String> parameters)
			throws IOException {
		super(parameters);
	}

	@Override
	public CompJob createCompJob(GenericJobMessage message, ToolboxTool tool, ResultCallback resultHandler) throws CompException {

		ToolDescription description = createToolDescription(tool);
		
		RCompJob analysisJob = new RCompJob();
		analysisJob.construct(message, description, resultHandler);
		analysisJob.setProcessPool(this.processPool);
		return analysisJob;
	}

	@Override
	protected String getStringDelimeter() {
		return RCompJob.STRING_DELIMETER;
	}

	@Override
	protected String getVariableNameSeparator() {
		return ".";
	}

}
