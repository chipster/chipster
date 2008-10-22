package fi.csc.microarray.client.workflow.api;

import fi.csc.microarray.client.operation.OperationDefinition;

public class WfOperationDefinition {

	private OperationDefinition operationDefinition;

	public WfOperationDefinition(OperationDefinition operationDefinition) {
		this.operationDefinition = operationDefinition;
	}

	public OperationDefinition getWrapped() {
		return operationDefinition;
	}

}
