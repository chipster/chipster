package fi.csc.microarray.analyser.java;

import fi.csc.microarray.MicroarrayException;

public class FailingAnalysisJob extends JavaAnalysisJobBase {

	@Override
	protected void execute() throws MicroarrayException {
		throw new MicroarrayException("Job failed.");
	}


	@Override
	public String getVVSADL() {		
		return " ANALYSIS Test/FailJava (Java job which fails.) ";
	}

}
