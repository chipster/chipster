package fi.csc.microarray.client.operation;

import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.client.operation.OperationRecord.ParameterRecord;
import fi.csc.microarray.client.operation.parameter.EnumParameter;
import fi.csc.microarray.client.operation.parameter.EnumParameter.SelectionOption;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.LinkUtils;
import fi.csc.microarray.description.GenericInputTypes;
import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLSyntax;
import fi.csc.microarray.description.SADLSyntax.InputType;
import fi.csc.microarray.module.chipster.ChipsterInputTypes;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.util.Strings;

/**
 * This class represents the "operations" that an user can select from the right
 * side list in the ToolSelectorPanel. These are the "blueprints" of specific
 * operations - the actual Operations
 * 
 * @author Janne Käki, Aleksi Kallio
 * 
 */
public class OperationDefinition implements ExecutionItem {

	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(OperationDefinition.class);

	/**
	 * Used to create an operation of origin for imported raw data. Not an
	 * actual operation (in the sense that it would be possible to execute it)
	 * but rather a dummy substitute, without any parameters.
	 */
	public static final String IMPORT_DEFINITION_ID = "operation-definition-id-import";
	public static final String CREATE_DEFINITION_ID = "operation-definition-id-user-modification";
	
	public static final OperationDefinition IMPORT_DEFINITION;
	public static final OperationDefinition CREATE_DEFINITION;

	/**
	 * An enumeration containing all possible results when evaluating an
	 * operation's suitability to a dataset.
	 */
	public static enum Suitability {
		SUITABLE, IMPOSSIBLE, EMPTY_REQUIRED_PARAMETERS;

		public boolean isOk() {
			return this == SUITABLE;
		}
	};

	/**
	 * The result of attempted binding. If binding was not successfull, then
	 * bindings is null.
	 */
	public static class BindingResult {
		Suitability suitability;
		LinkedList<DataBinding> bindings = null;
	}
	
	public static String IDENTIFIER_SEPARATOR = "/";

	
	static {
		// done here to guarantee right execution order
		IMPORT_DEFINITION = new OperationDefinition(IMPORT_DEFINITION_ID, "Import data",
		        ToolCategory.IMPORT_CATEGORY, "Import data.",
		        false, null);
		CREATE_DEFINITION = new OperationDefinition(CREATE_DEFINITION_ID, "Create a dataset",
		        ToolCategory.CREATE_CATEGORY, "Create a new dataset.",
		        false, null);
	}

	public static class InputDefinition {

		private String id;
		private String displayName;
		private String description = null;
		private String postfix = null;
		private boolean isMulti = false;
		private int multiCounter;
		private SADLSyntax.InputType type;
		private boolean isOptional;

		/**
		 * Creates single input.
		 * @param isOptional 
		 */
		public InputDefinition(Name name, String description, SADLSyntax.InputType type, boolean isOptional) {
			resetMulti();
			this.id = name.getID();
			this.displayName = name.getDisplayName();
			this.description = description;
			this.type = type;
			this.isOptional = isOptional;
		}

		/**
		 * Creates isMulti-input.
		 * @param isOptional 
		 */
		public InputDefinition(String prefix, String postfix, String displayName, String description, SADLSyntax.InputType type, boolean isOptional) {
			this.id = prefix;
			this.postfix = postfix;
			this.displayName = displayName;
			this.description = description;
			this.type = type;
			this.isOptional = isOptional;
			this.isMulti = true;
		}

		public String getID() {
			if (!isMulti) {
				return id;
			} else {
				return id + Strings.toString(multiCounter, 3) + postfix; // show always at least 3 digits 
			}
		}
		
        public String getDescription() {
                return this.description;
        }

        public String getDisplayName() {
        	return this.displayName;
        }
        
        /**
         * Get display name for isMulti input.
         * 
         * @param id
         * @return display name with {...} replaced with the number of the parameter id
         */
        public String getDisplayName(String id) {
        	if (!this.isMulti()) {
        		return this.displayName;
        	}
        	
        	if (!this.idMatches(id)) {
        		return this.displayName;
        	} 
        	
        	String middle = id.substring(this.id.length(), id.lastIndexOf(this.postfix));
        	return this.displayName.replaceAll("\\{\\.\\.\\.\\}", middle);
        }
        
        public void setDescription(String description) {
            this.description = description;
        }
		
		public SADLSyntax.InputType getType() {
		    return type;
		}

		public boolean isOptional() {
			return isOptional;
		}

		private void nextMulti() {
			multiCounter++;
		}

		public boolean isMulti() {
			return isMulti;
		}

		public void resetMulti() {
			multiCounter = 1;
		}
		
		public boolean idMatches(String id) {
			if (!this.isMulti()) {
				return this.getID().equals(id);
			} else {
				if (id.startsWith(this.id) && id.endsWith(this.postfix)) {
					String middle = id.substring(this.id.length(), id.lastIndexOf(this.postfix));
					if (middle.matches("\\d+")) {
						return true;
					}
				}
				return false;
			}
		}
	}

	private String id;
	private String displayName;
	private ToolCategory category;
	private LinkedList<Parameter> parameters = new LinkedList<Parameter>();
	private String description;
	private String helpURL;
	private int colorCount;
	private int outputCount = 0;
	private LinkedList<InputDefinition> inputs = new LinkedList<InputDefinition>();

	private boolean hasSourceCode;

	private boolean isLocal = false;

	/**
	 * Creates a new operation definition with the given initial values.
	 * 
	 * @param id
	 *            The name of this operation. Should be something that extends
	 *            the corresponding category name to be more specific (for
	 *            example, in the category "Normalization", "Lowess" might be a
	 *            good name of an operation definition, resulting in "Lowess
	 *            Normalization" when an actual operation is created).
	 * @param description
	 *            A written description of this operation's purpose.
	 */
	public OperationDefinition(String id, String displayName, ToolCategory category,
	                           String description, boolean hasSourceCode,
	                           String helpURL) {
		this.id = id;
		this.displayName = displayName;
		this.category = category;
		this.hasSourceCode = hasSourceCode;
		this.helpURL = helpURL;
		if (category != null) {
			category.addOperation(this);
		}

		this.description = description;
	}
	
	/**
	 * Like the constructor above, but makes it possible to define operation definition as local.
	 * 
	 * @param id
	 * @param displayName
	 * @param category
	 * @param description
	 * @param hasSourceCode
	 * @param helpURL
	 * @param isLocal
	 */
	public OperationDefinition(String id, String displayName, ToolCategory category,
            String description, boolean hasSourceCode,
            String helpURL, boolean isLocal) {
		this(id, displayName, category, description, hasSourceCode, helpURL);
		this.isLocal = isLocal;
	}

	/**
	 * Simplified constructor.
	 * @param id
	 * @param category
	 * @param description
	 * @param hasSourceCode
	 * @param helpURL
	 */
     public OperationDefinition(String id, String displayName, ToolCategory category,
         String description, boolean hasSourceCode) {
         this(id, displayName, category, description, hasSourceCode, null);
     }

	/**
	 * @return The name of this operation definition.
	 */
	public String getID() {
		return id;
	}

	public String getDisplayName() {
		return displayName;
	}
	
	/**
	 * @return The category to which this operation definition belongs.
	 */
	public ToolCategory getCategory() {
		return category;
	}

	/**
	 * @return The name of the category to which this belongs.
	 */
	public String getCategoryName() {
		return category.getName();
	}

	/**
	 * 
	 * @return categoryName / operationName
	 */
	public String getFullName() {
		return getCategoryName() + " / " + getDisplayName();
	}
	
	/**
     * @return URL linking to a help page or null if not given.
     */
    public String getHelpURL() {
        return helpURL;
    }

	/**
	 * @return A written description of this operation's purpose and function.
	 */
	public String getDescription() {
		return description;
	}

	/**
	 * @return A String representation (actually, just the name) of this
	 *         operation definition, used for showing this on the list.
	 */
	public String toString() {
		return getDisplayName();
	}

	/**
	 * Evaluates the suitability of this operation for the given dataset.
	 * 
	 * @param data
	 *            The dataset for which to evaluate.
	 * @param  
	 * @return One of the OperationDefinition.Suitability enumeration, depending
	 *         on how suitable the operation is judged.
	 */
	@Override
	public Suitability evaluateSuitabilityFor(Iterable<DataBean> data, List<DataBinding> bindings) {
	       
		if (bindings == null) {
			// Attempt binding
			BindingResult result = bindInputs(data);
			return result.suitability;
		} else {
			return evaluateBindingSuitability(data, bindings);
		}	    
	}
	
	public Suitability evaluateBindingSuitability(Iterable<DataBean> datas, List<DataBinding> bindings) {
		
		// check that all mandatory inputs are set
		for (InputDefinition input : inputs) {
			if (!input.isOptional()) {
				if (!bindingsContainsInput(bindings, input)) {
					logger.debug("  no binding found for " + input.getID());
					return Suitability.EMPTY_REQUIRED_PARAMETERS;
				}
			}
		}
		
		// collect a set of bound datasets
		HashSet<DataBean> boundDatas = new HashSet<>();  
		for (DataBinding binding : bindings) {
			boundDatas.add(binding.getData());
		}	

		/*
		 * Special case: Don't care about selected datasets if the tools doesn't
		 * take any inputs, because in this case it's not obvious that users
		 * should clear the selection.
		 */
		if (!inputs.isEmpty()) {
			// check that all selected datasets are bound
			for (DataBean data : datas) {
				if (!boundDatas.contains(data)) {
					logger.debug("  concrete input " + data.getName() + " was not bound");
					return Suitability.EMPTY_REQUIRED_PARAMETERS;
				}
			}
		}
        
        return Suitability.SUITABLE;
	}
	
	private boolean bindingsContainsInput(List<DataBinding> bindings,
			InputDefinition input) {

		for (DataBinding binding : bindings) {
			if (input.idMatches(binding.getName())) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Check suitability of a given parameter list. The parameter
	 * list can also come from the Operation object that encapsulates
	 * this definition.
	 * 
	 * @param params
	 * @return
	 */
	public static Suitability evaluateParameterSuitability(List<Parameter> params) {
        for (Parameter param : params) {
            // Required parameters can not be empty
            if (!param.isOptional() && (param.getValue() == null ||
                                        param.getValue().equals(""))) {
                return Suitability.EMPTY_REQUIRED_PARAMETERS;
            }
        }
        
        return Suitability.SUITABLE;
	}

	public LinkedList<Parameter> getParameters() {
		return parameters;
	}
	
	public Parameter getParameter(String id) {
		for (Parameter parameter: parameters) {
			if (parameter.getID().equals(id)) {
				return parameter;
			}
		}
		return null;
	}

	public void addParameter(Parameter parameter) {
		parameters.add(parameter);
	}

	public int getColorCount() {
		return colorCount;
	}

	public void addInput(Name name, String description, InputType type, boolean isOptional) {
		InputDefinition input = new InputDefinition(name, description, type, isOptional);
		inputs.add(input);
	}

	public void addInput(String prefix, String postfix, String displayName, String description, InputType type, boolean isOptional) {
		InputDefinition input = new InputDefinition(prefix, postfix, displayName, description, type, isOptional);
		inputs.add(input);
	}
	
	public List<InputDefinition> getInputs() {
	    return inputs;
	}

	/**
	 * 
	 * @param id
	 * @return null if no input with the given id is found
	 */
	public InputDefinition getInput(String id) {
		for (InputDefinition input : inputs) {
			if (input.idMatches(id)) {
				return input;
			}
		}
		return null;
	}
	
	/**
	 * In a nutshell, formal inputs (as defined by the operation) are bound to
	 * concrete inputs (as chosen by user) using greedy and order-based
	 * algorithm. Formal inputs are processed in order and first fitting
	 * concrete input is bound to them. If formal input can have multiple
	 * concrete inputs, then all fitting ones are bound. Always at least one
	 * concrete input must be bound, a single concrete input cannot be bound
	 * multiple times and in the end all concrete inputs must be bound.
	 * 
	 * @param inputValues
	 *            no changes are made to this parameter
	 * @return BindingResult object, where bindings is null when binding failed
	 */
	public BindingResult bindInputs(Iterable<DataBean> inputValues) {

		BindingResult result = new BindingResult();
		
		// initialise
		LinkedList<DataBinding> bindings = new LinkedList<DataBinding>();
		LinkedList<DataBean> notProcessedInputValues = new LinkedList<DataBean>();
		for (DataBean bean : inputValues) {
			notProcessedInputValues.add(bean);
		}

		LinkedList<InputDefinition> unboundMetadataDefinitions = new LinkedList<InputDefinition>();

		logger.debug("binding " + notProcessedInputValues.size() + " values to " + inputs.size() + " formal inputs");

		// bind by iterating over formal parameters
		for (InputDefinition input : inputs) {
			input.resetMulti();

			// metadata needs not to be selected, it is fetched automatically
			if (doBackwardsCompatibleMetadataCheck(input)) {
				unboundMetadataDefinitions.add(input);
				continue;
				
			}

			// find values to bind by iterating over remaining actual parameters
			LinkedList<DataBean> removedValues = new LinkedList<DataBean>();
			for (DataBean value : notProcessedInputValues) {

				// try to match values to input definitions
				logger.debug("  trying to bind " + value.getName() + " to " + input.id + " (" + input.type + ")");
				if (doBackwardsCompatibleTypeCheck(input.type, value)) {

					logger.debug("    bound successfully (" + value.getName() + " -> " + input.getID() + ")");

					bindings.add(new DataBinding(value, input.getID(), input.getType()));
					removedValues.add(value); // mark it to be removed after iteration
					
					if (!input.isMulti()) {
						break;
					} else {
						input.nextMulti();
					}
				}
			}
			notProcessedInputValues.removeAll(removedValues);
		}
		

		// automatically bind phenodata, if needed
		logger.debug("we have " + bindings.size() + " bindings before metadata retrieval");
		if (!unboundMetadataDefinitions.isEmpty()) {

			LinkedList<DataBinding> phenodataBindings = new LinkedList<DataBinding>(); // need this to prevent ConcurrentModificationException
			
			Iterator<DataBinding> bindingIterator = bindings.iterator();
			Iterator<InputDefinition> unboundMetadataIterator= unboundMetadataDefinitions.iterator();
						
			while (bindingIterator.hasNext() && unboundMetadataIterator.hasNext()) {
				
				DataBean input = bindingIterator.next().getData(); // bind inputs and phenodatas in same order
				InputDefinition unboundMetadata = unboundMetadataIterator.next();
				
				// locate annotation (metadata) link from input bean or one of its ancestors				
				DataBean metadata = LinkUtils.retrieveInherited(input, Link.ANNOTATION);

				if (metadata != null) {
					phenodataBindings.add(new DataBinding(metadata, unboundMetadata.getID(), ChipsterInputTypes.PHENODATA));
					
				} else {
					result.suitability = Suitability.IMPOSSIBLE;
					return result;
				}
			}
			bindings.addAll(phenodataBindings);
		}
		
		logger.debug("we have " + bindings.size() + " bindings after metadata retrieval");
		
		if (!evaluateBindingSuitability(inputValues, bindings).isOk()) {
			result.suitability = Suitability.IMPOSSIBLE;
			return result;
		}

		// return successful binding result
		result.bindings = bindings;
		result.suitability = Suitability.SUITABLE;
		return result;
	}

	// TODO update to new type tag system
	private boolean doBackwardsCompatibleTypeCheck(InputType type, DataBean data) {
		
		if (type == ChipsterInputTypes.AFFY) {
			return data.hasTypeTag(MicroarrayModule.TypeTags.RAW_AFFYMETRIX_EXPRESSION_VALUES);
			
		} else if (type == ChipsterInputTypes.CDNA) {
			return data.hasTypeTag(MicroarrayModule.TypeTags.RAW_EXPRESSION_VALUES);
			
		} else if (type == ChipsterInputTypes.GENE_EXPRS) {
			return data.hasTypeTag(MicroarrayModule.TypeTags.NORMALISED_EXPRESSION_VALUES);
			
		} else if (type == ChipsterInputTypes.GENELIST) {
			return data.hasTypeTag(MicroarrayModule.TypeTags.GENENAMES);
		
		} else if (type == ChipsterInputTypes.BAM) {
			return data.hasTypeTag(MicroarrayModule.TypeTags.BAM_FILE);
			
		} else if (type == ChipsterInputTypes.FASTA) {
			return data.hasTypeTag(MicroarrayModule.TypeTags.FASTA_FILE);
			
		} else if (type == ChipsterInputTypes.GTF) {
			return data.hasTypeTag(MicroarrayModule.TypeTags.GTF_FILE);
			
		} else if (type == ChipsterInputTypes.PHENODATA) {
			return data.hasTypeTag(MicroarrayModule.TypeTags.PHENODATA);
			
		} else if (type == GenericInputTypes.GENERIC) {
			return true;
			
		} else if (type == ChipsterInputTypes.MOTHUR_OLIGOS) {
			return data.hasTypeTag(MicroarrayModule.TypeTags.MOTHUR_OLIGOS);
			
		} else if (type == ChipsterInputTypes.MOTHUR_NAMES) {
			return data.hasTypeTag(MicroarrayModule.TypeTags.MOTHUR_NAMES);
			
		} else if (type == ChipsterInputTypes.MOTHUR_GROUPS) {
				return data.hasTypeTag(MicroarrayModule.TypeTags.MOTHUR_GROUPS);
		
		} else {
			throw new IllegalArgumentException();
		}
	}

	// TODO update to meta input type
	private boolean doBackwardsCompatibleMetadataCheck(InputDefinition input) {
		return input.id.startsWith("phenodata");
	}

	public int getOutputCount() {
		return this.outputCount;
	}

	public void setOutputCount(int outputCount) {
		this.outputCount = outputCount;
	}

	public boolean hasSourceCode() {
		return hasSourceCode;
	}
	
	public String toStringVerbose() {
		String s = "\n-------------- operation definition --------------\n";
		s += getCategoryName() + " / ";
		s += getDisplayName() + " ";
		s += "(" + getID() + ")\n";
		for (InputDefinition input: inputs) {
			String type;
			if (input.getType() != null) {
				type = input.getType().getName();
			} else {
				type = "null";
			}
			s += input.getID() + " " + type + " " + input.getDescription() + "\n";
		}
		for (Parameter parameter: parameters) {
		    // Some parameters don't have default values
		    String value;
		    if (parameter.getValue() == null) {
		        value = "[no default value]";
		    } else {
		        value = parameter.getValueAsString();
		    }
			s += parameter.getID() + " " + value + "\n";
		}

		s += "\n-------------- operation definition --------------\n";
		
		return s;
	}

	/**
	 * @return true if this operation should be run in LocalTaskExecutor instead of server
	 */
	public boolean isLocal() {
		return isLocal;
	}
	
	/**
	 * Checks if the operation defined by this operation definition is safe to run as batch.
	 * Batch run is tricky for input binding and input sensitive parameter evaluation, so 
	 * we do strict checks and only allow "simple" operations to be run as batch. This is 
	 * done to not confuse the user and create surprising effects to workflows, for example.
	 * 
	 * @return true if the operation passes checks and is safe to run as batch
	 */
	public boolean isBatchable() {
		
		// Check that parameters are not sensitive to inputs
		for (Parameter parameter : parameters) {
			if (parameter.isInputSensitive()) {
				return false;
			}
		}
		
		// Allow only single input operations to be batched
		return inputs.size() == 1;
	}

	public String getParameterDefaultValue(ParameterRecord parameterRecord) {
		Parameter parameter = getParameter(parameterRecord.getNameID().getID());					
		String defaultValue = null;
		if (parameter != null) {
			 defaultValue  = parameter.getValueAsString();
		}
		return defaultValue;
	}

	public String getHumanReadableParameterValue(ParameterRecord parameterRecord) {
		Parameter parameter = getParameter(parameterRecord.getNameID().getID());
		
		//EnumParameters have a display name for values
		if (parameter instanceof EnumParameter) {
			EnumParameter enumParameter = (EnumParameter) parameter;
			
			Object[] options = enumParameter.getOptions();
			
			// column selection doesn't have better name
			if (options != null) {
				// search for human readable name
				for (Object choice : options) {
					SelectionOption option = (SelectionOption) choice;
					
					//for (SelectionOption option : options) {
					if (parameterRecord.getValue().equals(option.getValue())) {
						
						return option.toString();
					}
				}
			}
		}
		return parameterRecord.getValue();
	}
}
