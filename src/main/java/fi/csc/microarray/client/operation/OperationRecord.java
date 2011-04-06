package fi.csc.microarray.client.operation;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.databeans.DataBean;

public class OperationRecord {

	private String id;
	private String displayName;
	
	private String sourceCode;

	private LinkedList<Parameter> parameters;
	private HashMap<String, DataBean> inputs;
	
	public OperationRecord(Operation operation) {
		
		// name
		this.id = operation.getID();
		this.displayName = operation.getDisplayName();
	
		// parameters
		for (Parameter parameter : operation.getParameters()) {
			this.parameters.add((Parameter) parameter.clone());
		}
	
		// inputs
		for (DataBinding binding : operation.getBindings()) {
			this.inputs.put(binding.getName(), binding.getData());
		}
	}
	
	public String getID() {
		return id;
	}
	
	public String getDisplayName() {
		return displayName;
	}

	public List<Parameter> getParameters() {
		return parameters;
	}
	
	
}
