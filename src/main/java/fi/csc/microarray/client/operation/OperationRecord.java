package fi.csc.microarray.client.operation;

import java.awt.Color;
import java.util.Collection;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.LinkedList;

import fi.csc.microarray.client.NameID;
import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.client.operation.OperationDefinition.InputDefinition;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.databeans.DataBean;


/**
 * Stores information about an executed Operation.
 * 
 * 
 * @author hupponen
 *
 */
public class OperationRecord {

	private NameID nameID = new NameID();;

	private String categoryName;
	private Color categoryColor;
	
	private String moduleName;
	
	private String sourceCode;

	private LinkedHashMap<String, ParameterRecord> parameters = new LinkedHashMap<String, ParameterRecord>();
	private LinkedHashMap<String, InputRecord> inputs = new LinkedHashMap<String, InputRecord>();

	// jobId is set for operations that are still running and null for completed
	// operations
	private String jobId;
	
	private Date startTime;
	private Date endTime;
	
	public Date getStartTime() {
		return startTime;
	}

	public void setStartTime(Date startTime) {
		this.startTime = startTime;
	}

	public Date getEndTime() {
		return endTime;
	}

	public void setEndTime(Date endTime) {
		this.endTime = endTime;
	}

	public OperationRecord(Operation operation) {
		
		// name
		this.nameID.setId(operation.getID());
		this.nameID.setDisplayName(operation.getDisplayName());
		this.nameID.setDescription(operation.getDescription());
	
		// category
		this.categoryName = operation.getCategoryName();
		this.categoryColor = operation.getCategoryColor();
		
		// module
		if (operation.getDefinition().getCategory().getModule() != null) {
			this.moduleName = operation.getDefinition().getCategory().getModule().getModuleName();
		}
		
		// parameters
		for (Parameter parameter : operation.getParameters()) {
			this.parameters.put(parameter.getID(), new ParameterRecord(parameter));
		}
	
		// inputs
		for (DataBinding binding : operation.getBindings()) {
			InputDefinition inputDefinition = operation.getDefinition().getInput(binding.getName());
			String displayName = inputDefinition.isMulti() ? inputDefinition.getDisplayName(binding.getName()) : inputDefinition.getDisplayName();
			this.inputs.put(binding.getName(), new InputRecord(binding.getName(), displayName, inputDefinition.getDescription(), binding.getData()));
		}
	
		// source code
		// TODO get it from ResultMessage when it becomes available
		this.sourceCode = "not yet available";
	
	}
	
	public OperationRecord() {
	}

	public NameID getNameID() {
		return this.nameID;
	}
	
	public String getCategoryName() {
		return categoryName;
	}

	public String getModule() {
		return this.moduleName;
	}
	
	public String getFullName() {
		return getCategoryName() + " / " + nameID.getDisplayName();
	}

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

	public ParameterRecord getParameter(String id) {
		return parameters.get(id);
	}

	/**
	 * Iteration order of the parameters is determined.
	 * @return
	 */
	public Collection<InputRecord> getInputRecords() {
		return inputs.values();
	}

	public String getSourceCode() {
		return sourceCode;
	}

	public void setCategoryName(String categoryName) {
		this.categoryName = categoryName;
	}

	public void setModule(String moduleName) {
		this.moduleName = moduleName;
	}
	
	public void setCategoryColor(Color categoryColor) {
		this.categoryColor = categoryColor;
	}

	public void setSourceCode(String sourceCode) {
		this.sourceCode = sourceCode;
	}


	public void addParameter(Parameter parameter) {
		this.parameters.put(parameter.getID(), new ParameterRecord(parameter));
	}
	
	public void addParameter(String id, String displayName, String description, String value) {
		this.parameters.put(id, new ParameterRecord(id, displayName, description, value));
	}

	
	public void addInput(NameID nameID, DataBean dataBean) {
		this.inputs.put(nameID.getID(), new InputRecord(nameID, dataBean));
	}
	
	public void addInput(NameID nameId, String dataId) {
		this.inputs.put(nameID.getID(), new InputRecord(nameID, dataId));
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


	private class StringRecord {
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

	private class DataBeanRecord {
		private NameID nameID;
		private String dataId; // for cases where DataBean itself is missing
		private DataBean value;

		public DataBeanRecord(NameID nameID, DataBean value) {
			this.nameID = nameID;
			if (value != null) {
				dataId = value.getId();
			}
			this.value = value;
		}
		
		public DataBeanRecord(String id, String displayName, String description, DataBean value) {
			this(new NameID(id, displayName, description), value);
		}
		
		public DataBeanRecord(NameID nameID, String dataId) {
			this.nameID = nameID;
			this.dataId = dataId;
			this.value = null;
		}
		
		
		public NameID getNameID() {
			return nameID;
		}

		public String getDataId() {
			return dataId;
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
		
		public InputRecord (NameID nameID, String dataId) {
			super(nameID, dataId);
		}
		
	}

	public void setNameID(NameID nameID) {
		this.nameID = nameID;
	}

	public static OperationRecord getUnkownOperationRecord() {
		OperationRecord record = new OperationRecord();
		record.setNameID(new NameID("unknown", "Unknown", "No tool information available."));
		record.setCategoryName("Unknown");
		record.setCategoryColor(ToolCategory.UNKNOWN_CATEGORY_COLOR);
		record.setSourceCode("Source code not available.");
		return record;
	}

	public Iterable<DataBean> getInputDataBeans() {
		LinkedList<DataBean> beans = new LinkedList<DataBean>();
		for (InputRecord input : getInputRecords()) {
			beans.add(input.getValue());
		}
		return beans;
	}

	public void setJobId(String jobId) {
		this.jobId = jobId;
	}
	
	public String getJobId() {
		return this.jobId;
	}
}
