package fi.csc.microarray.comp.java;

import fi.csc.microarray.comp.JobCancelledException;
import fi.csc.microarray.messaging.JobState;

public class FailingCompJob extends JavaCompJobBase {

	@Override
	protected void execute() throws JobCancelledException {
		outputMessage.setErrorMessage("This job always fails.");
		outputMessage.setOutputText("There's no way around this.");
		updateState(JobState.FAILED, "");
	}


	@Override
	public String getSADL() {		
		return " ANALYSIS Test/FailJava (Java job which fails.) ";
	}

}
