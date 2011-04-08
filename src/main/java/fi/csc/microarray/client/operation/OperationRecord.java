package fi.csc.microarray.client.operation;

import java.awt.Color;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;

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
	
		// source code
		this.sourceCode = "not yet available";
	
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
	
	public String getSourceCode() {
		return sourceCode;
	}

	public Parameter getParameter(String name) {
		for (Parameter parameter : parameters) {
			if (parameter.getID().equals(name)) {
				return parameter;
			}
		}
		return null;
	}

	/**
	 * If the given DataBean is a value of an input binding in this record,
	 * set the value of such input binding to null.
	 * 
	 * @param inputBean
	 */
	public void removeInput(DataBean inputBean) {
		if (inputBean == null) {
			return;
		}
		for (Entry<String, DataBean> entry : inputs.entrySet()) {
			if (entry.getValue() == inputBean) {
				entry.setValue(null);
			}
		}
	}

}
