package fi.csc.microarray.client.workflow.api;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.operation.OperationDefinition;

public class WfApplication {

	private ClientApplication application;

	public WfApplication(ClientApplication application) {
		this.application = application;
	}

	public WfSelectionManager getSelectionManager() {
		return new WfSelectionManager(application.getSelectionManager());
	}
	
	public void executeOperation(final WfOperation operation) {
		this.application.executeOperation(operation.getWrapped());
	}
	
	public WfOperationDefinition getOperationDefinition(String operationID) {
		OperationDefinition operationDefinition = this.application.getOperationDefinition(operationID);
		return new WfOperationDefinition(operationDefinition);
	}

}
