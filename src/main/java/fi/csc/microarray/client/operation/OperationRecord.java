package fi.csc.microarray.client.operation;

import java.awt.Color;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.databeans.DataBean;

public class OperationRecord {

	private String id;
	private String displayName;
	
	private String categoryName;
	private Color categoryColor;
	
	private String sourceCode;

	private List<Parameter> parameters = new LinkedList<Parameter>();
	private LinkedHashMap<String, DataBean> inputs = new LinkedHashMap<String, DataBean>();
	
	
	public OperationRecord(Operation operation) {
		
		// name
		this.id = operation.getID();
		this.displayName = operation.getDisplayName();
	
		// category
		this.categoryName = operation.getCategoryName();
		this.categoryColor = operation.getCategoryColor();
		
		// parameters
		for (Parameter parameter : operation.getParameters()) {
			this.parameters.add((Parameter) parameter.clone());
		}
	
		// inputs
		if (operation.getBindings() != null) {
			for (DataBinding binding : operation.getBindings()) {
				this.inputs.put(binding.getName(), binding.getData());
			}
		}
	}
	
	public String getID() {
		return id;
	}
	
	public String getDisplayName() {
		if (displayName != null && !displayName.isEmpty()) {
			return displayName;
		} else {
			return getID();
		}
	}

	public String getCategoryName() {
		return categoryName;
	}
	
	public String getFullName() {
		return getCategoryName() + " / " + getDisplayName();
	}
	
	
	// FIXME check if operation with the id still exists and if the color has changed
	public Color getCategoryColor() {
		return categoryColor;
	}
	
	public List<Parameter> getParameters() {
		return parameters;
	}
	
	public LinkedHashMap<String, DataBean> getInputs() {
		return inputs;
	}
	
	
}
