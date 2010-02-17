package fi.csc.microarray.analyser.emboss;

import fi.csc.microarray.analyser.AnalysisDescription;
import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.messaging.message.ResultMessage;

public class EmbossAnalysisJob extends AnalysisJob {

	@Override
	protected void cancelRequested() {
		// TODO Auto-generated method stub

	}

	@Override
	protected void execute() throws Exception {
		// TODO this method is just a stub for now
		
		// this is our input:
		JobMessage jobMessage = this.inputMessage;
		
		// this is our analysis description (compute service specific version 
		// of SADLDescription):
		AnalysisDescription description = this.analysis;
		
		// this is what we should produce as output:
		ResultMessage outputMessage = this.outputMessage;
		
		// and this is where results are returned 
		ResultCallback resultHandler = this.resultHandler;
		
		// there are our parameters (their order is significant: nth value corresponds 
		// to nth parameter in AnalysisDescription):
		for (String paramValue : this.inputMessage.getParameters()) {
			System.out.println(paramValue);
		}

	}

}
