package fi.csc.microarray.analyser.java;

import fi.csc.microarray.analyser.JobCancelledException;
import fi.csc.microarray.messaging.JobState;

public class FailingAnalysisJob extends JavaAnalysisJobBase {

	@Override
	protected void execute() throws JobCancelledException {
		outputMessage.setErrorMessage("This job always fails.");
		outputMessage.setOutputText("There's no way around this.");
		updateState(JobState.FAILED, "", true);
	}


	@Override
	public String getVVSADL() {		
		return " ANALYSIS Test/FailJava (Java job which fails.) ";
	}

}
