package fi.csc.microarray.analyser.java;

import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.JobCancelledException;

public abstract class JavaAnalysisJobBase extends AnalysisJob {

	@Override
	protected void preExecute() throws JobCancelledException  {
		super.preExecute();		
	}

	@Override
	protected void postExecute() throws Exception {
		super.postExecute();
	}

	@Override
	protected void cleanUp() {
		super.cleanUp();
	}
	
	@Override
	protected void cancelRequested() {
		// ignore by default
	}

	public abstract String getVVSADL();
}
