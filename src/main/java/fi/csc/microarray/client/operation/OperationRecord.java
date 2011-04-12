package fi.csc.microarray.client.operation;

import java.awt.Color;
import java.util.Collection;
import java.util.LinkedHashMap;

import fi.csc.microarray.client.NameID;
import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.client.operation.OperationDefinition.InputDefinition;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.databeans.DataBean;

public class OperationRecord {

	
	private String id;
	private String displayName;
	private String description;
	
	public String getDescription() {
		return description;
	}

	public void setDescription(String description) {
		this.description = description;
	}

	private String categoryName;
	private Color categoryColor;
	
	private String sourceCode;

	private LinkedHashMap<String, ParameterRecord> parameters = new LinkedHashMap<String, ParameterRecord>();
	private LinkedHashMap<String, InputRecord> inputs = new LinkedHashMap<String, InputRecord>();
	
	
	public OperationRecord(Operation operation) {
		
		// name
		this.id = operation.getID();
		this.displayName = operation.getDisplayName();
	
		// category
		this.categoryName = operation.getCategoryName();
		this.categoryColor = operation.getCategoryColor();
		
		// parameters
		for (Parameter parameter : operation.getParameters()) {
			this.parameters.put(parameter.getID(), new ParameterRecord(parameter));
		}
	
		// inputs
		if (operation.getBindings() != null) {
			for (DataBinding binding : operation.getBindings()) {
				InputDefinition inputDefinition = operation.getDefinition().getInput(binding.getName());
				String displayName = inputDefinition.isMulti() ? inputDefinition.getDisplayName(binding.getName()) : inputDefinition.getDisplayName();
				this.inputs.put(binding.getName(), new InputRecord(binding.getName(), displayName, inputDefinition.getDescription(), binding.getData()));
			}
		}
	
		// source code
		this.sourceCode = "not yet available";
	
	}
	
	public OperationRecord() {
	}

	public void setID(String id) {
		this.id = id;
	}

	public void setDisplayName(String displayName) {
		this.displayName = displayName;
	}

	public void setCategoryName(String categoryName) {
		this.categoryName = categoryName;
	}

	public void setCategoryColor(Color categoryColor) {
		this.categoryColor = categoryColor;
	}

	public void setSourceCode(String sourceCode) {
		this.sourceCode = sourceCode;
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
	
	/**
	 * Iteration order of the parameters is determined.
	 * @return
	 */
	public Collection<ParameterRecord> getParameters() {
		return parameters.values();
	}

	/**
	 * Iteration order of the parameters is determined.
	 * @return
	 */
	public Collection<InputRecord> getInputs() {
		return inputs.values();
	}
	
	public String getSourceCode() {
		return sourceCode;
	}

	public ParameterRecord getParameter(String id) {
		return parameters.get(id);
	}

	public void addParameter(Parameter parameter) {
		this.parameters.put(parameter.getID(), new ParameterRecord(parameter));
	}
	
	public void addParameter(String id, String displayName, String description, String value) {
		this.parameters.put(id, new ParameterRecord(id, displayName, description, value));
	}

	
	public void addInput(NameID nameID, DataBean dataBean) {
		this.inputs.put(id, new InputRecord(nameID, dataBean));
	}
	
	
	/**
	 * If the given DataBean is a value of an InputRecord in this OperationRecord,
	 * set the value of such input binding to null.
	 * 
	 * @param inputBean
	 */
	public void removeInput(DataBean inputBean) {
		if (inputBean == null) {
			return;
		}
		for (InputRecord inputRecord : inputs.values()) {
			if (inputRecord.getValue() == inputBean) {
				inputRecord.setValue(null);
			}
		}
	}


	public class StringRecord {
		private NameID nameID;
		private String value;

		public StringRecord(String id, String displayName, String description, String value) {
			this.nameID = new NameID(id, displayName, description);
			this.value = value;
		}
		
		public NameID getNameID() {
			return nameID;
		}

		public String getValue() {
			return value;
		}
	}

	public class ParameterRecord extends StringRecord {
		
		public ParameterRecord(String id, String displayName, String description, String value) {
			super(id, displayName, description, value);
		}
		
		public ParameterRecord(Parameter parameter) {
			super(parameter.getID(), parameter.getDisplayName(), parameter.getDescription(), parameter.getValueAsString());
		}
	}

	public class DataBeanRecord {
		private NameID nameID;
		private DataBean value;

		public DataBeanRecord(NameID nameID, DataBean value) {
			this.nameID = nameID;
			this.value = value;
		}
		
		public DataBeanRecord(String id, String displayName, String description, DataBean value) {
			this.nameID = new NameID(id, displayName, description);
			this.value = value;
		}
		
		public NameID getNameID() {
			return nameID;
		}

		public DataBean getValue() {
			return value;
		}

		public void setValue(DataBean value) {
			this.value = value;
		}

	}
	
	public class InputRecord extends DataBeanRecord {

		public InputRecord (NameID nameID, DataBean value) {
			super(nameID, value);
		}
		
		public InputRecord(String id, String displayName, String description, DataBean value) {
			super(id, displayName, description, value);
		}
	}
	
}
