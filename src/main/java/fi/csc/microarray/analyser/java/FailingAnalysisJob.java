package fi.csc.microarray.analyser.java;

import fi.csc.microarray.exception.MicroarrayException;

public class FailingAnalysisJob extends JavaAnalysisJobBase {

	@Override
	protected void execute() throws MicroarrayException {
		throw new MicroarrayException("Job failed.");
	}


	@Override
	public String getSADL() {		
		return " ANALYSIS Test/FailJava (Java job which fails.) ";
	}

}
