package fi.csc.microarray.client.operation;

import java.awt.Color;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.LinkUtils;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.description.GenericInputTypes;
import fi.csc.microarray.description.SADLSyntax;
import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLSyntax.InputType;
import fi.csc.microarray.module.basic.BasicModule;
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
	 * 
	 * @author Janne Käki
	 * 
	 */
	public static enum Suitability {
		SUITABLE, IMPOSSIBLE, ALREADY_DONE, TOO_MANY_INPUTS, NOT_ENOUGH_INPUTS,
		EMPTY_REQUIRED_PARAMETERS;

		private static final Color GREEN = new Color(52, 196, 49);
		private static final Color YELLOW = new Color(196, 186, 49);
		private static final Color RED = new Color(196, 49, 49);
		private static final Color COLOR_ALREADY_DONE = new Color(230, 180, 250);

		public boolean isImpossible() {
			return this == IMPOSSIBLE || this == NOT_ENOUGH_INPUTS || this == TOO_MANY_INPUTS;
		}

		public boolean isOk() {
			return this == SUITABLE;
		}

		/**
		 * @return The indicator color of this Suitability item (to be used, for
		 *         example, as the background of the suitability label).
		 */
		public Color getIndicatorColor() {
			if (isImpossible()) {
				return RED;
			} else if (isOk()) {
				return GREEN;
			} else if (this == ALREADY_DONE) {
				return COLOR_ALREADY_DONE;
			} else {
				return YELLOW;
			}
		}

		/**
		 * @return A String representation of this Suitability item (to be used,
		 *         for example, as the text in the suitability label).
		 */
		public String toString() {
			switch (this) {
			case SUITABLE:
				return "Suitable";
			case IMPOSSIBLE:
				return "Impossible";
			case ALREADY_DONE:
				return "Already done";
			case TOO_MANY_INPUTS:
				return "Too many inputs";
			case NOT_ENOUGH_INPUTS:
				return "Not enough inputs";
            case EMPTY_REQUIRED_PARAMETERS:
                return "Some required parameters are empty";
			default:
				throw new RuntimeException("unknown suitability: " + this.name());
			}
		}
	};

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
	private Suitability evaluatedSuitability = null;

	private boolean hasSourceCode;

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
	 * @param parametersSuitability is either null - indicating that the
	 *        parameter suitability has not been checked yet or Suitability
	 *        object defining the suitability of parameters in an encapsulating
	 *        Operation object that calls this method.
	 * @return One of the OperationDefinition.Suitability enumeration, depending
	 *         on how suitable the operation is judged.
	 */
	public Suitability evaluateSuitabilityFor(Iterable<DataBean> data,
	        Suitability parameterSuitability) {
	       
        // Input suitability gets checked while trying to bind the data
        bindInputs(data);
	    
	    // Report only either input or parameter suitability
	    if (evaluatedSuitability.isOk()) {
	        if (parameterSuitability == null) {
	        	evaluatedSuitability = parameterSuitability(getParameters());
	        } else {
	        	evaluatedSuitability = parameterSuitability;
	        }
	    }
		
		return getEvaluatedSuitability();
	}
	
	/**
	 * Check suitability of a given parameter list. The parameter
	 * list can also come from the Operation object that encapsulates
	 * this definition.
	 * 
	 * @param params
	 * @return
	 */
	public static Suitability parameterSuitability(List<Parameter> params) {
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
	 * @return null when binding could not be done
	 */
	public LinkedList<DataBinding> bindInputs(Iterable<DataBean> inputValues) {

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
			boolean foundBinding = false;

			// metadata needs not to be selected, it is fetched automatically
			if (doBackwardsCompatibleMetadataCheck(input)) {
				foundBinding = true; // we'll find it later on
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
					foundBinding = true;
					removedValues.add(value); // mark it to be removed after iteration
					
					if (!input.isMulti()) {
						break;
					} else {
						input.nextMulti();
					}
				}
			}
			notProcessedInputValues.removeAll(removedValues);

			// input not bound and is mandatory, so can give up
			if (!foundBinding && !input.isOptional()) {
				logger.debug("  no binding found for " + input.id);
				this.evaluatedSuitability = Suitability.NOT_ENOUGH_INPUTS;
				return null;
			}
		}
		if (notProcessedInputValues.size() > 0) {
			logger.debug("  " + notProcessedInputValues.size() + " concrete inputs were not bound");
			this.evaluatedSuitability = Suitability.TOO_MANY_INPUTS;
			return null;
		}

		// automatically bind phenodata, if needed
		logger.debug("we have " + bindings.size() + " bindings before metadata retrieval");
		if (!unboundMetadataDefinitions.isEmpty()) {

			Iterator<DataBinding> bindingIterator = bindings.iterator();
			LinkedList<DataBinding> phenodataBindings = new LinkedList<DataBinding>(); // need this to prevent ConcurrentModificationException
			for (InputDefinition unboundMetadata : unboundMetadataDefinitions) {
				
				// locate annotation (metadata) link from input bean or one of its ancestors				
				DataBean input = bindingIterator.next().getData(); // bind inputs and phenodatas in same order
				DataBean metadata = LinkUtils.retrieveInherited(input, Link.ANNOTATION);

				if (metadata != null) {
					phenodataBindings.add(new DataBinding(metadata, unboundMetadata.getID(), ChipsterInputTypes.PHENODATA));
					
				} else {
					this.evaluatedSuitability = Suitability.NOT_ENOUGH_INPUTS;
					return null;
				}
			}
			bindings.addAll(phenodataBindings);
		}		
		logger.debug("we have " + bindings.size() + " bindings after metadata retrieval");

		this.evaluatedSuitability = Suitability.SUITABLE;
		return bindings;
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
			return data.hasTypeTag(BasicModule.TypeTags.PHENODATA);
			
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

	// TODO update to new type tag system
	private boolean doBackwardsCompatibleMetadataCheck(InputDefinition input) {
		return input.id.startsWith("phenodata");
	}

	/**
	 * @return the suitability of last bindInputs-call or null
	 * @see #bindInputs(Iterable)
	 */
	public Suitability getEvaluatedSuitability() {
		return evaluatedSuitability;
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
	
}
