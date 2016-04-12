package fi.csc.microarray.comp.python;

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

/**
 * Handler for analysis tools written in Python.
 * 
 * @author Aleksi Kallio
 *
 */
public class PythonJobFactory extends InterpreterJobFactory {

	static final Logger logger = Logger
			.getLogger(PythonJobFactory.class);
	

	public PythonJobFactory(HashMap<String, String> parameters)
			throws IOException {
		super(parameters);
	}

	@Override
	public CompJob createCompJob(GenericJobMessage message, ToolboxTool tool, ResultCallback resultHandler) throws CompException {
		ToolDescription description = createToolDescription(tool);
		
		PythonCompJob analysisJob = new PythonCompJob();
		analysisJob.construct(message, description, resultHandler);
		analysisJob.setProcessPool(this.processPool);
		return analysisJob;
	}

	@Override
	protected String getStringDelimeter() {
		return PythonCompJob.STRING_DELIMETER;
	}


	@Override
	protected String getVariableNameSeparator() {
		return "_";
	}

}
