package fi.csc.microarray.description;

import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.description.SADLSyntax.InputType;
import fi.csc.microarray.description.SADLSyntax.ParameterType;

/**
 * SADL description for one analysis tool. Describes analysis tools so 
 * that they can be used in Chipster context.
 * 
 * @author Aleksi Kallio
 * 
 */
public class SADLDescription {

	private Name name;
	private String comment;
	
	private LinkedList<Input> inputs = new LinkedList<Input>();
	private LinkedList<Output> outputs = new LinkedList<Output>();
	private LinkedList<Input> metaInputs = new LinkedList<Input>();
	private LinkedList<Output> metaOutputs = new LinkedList<Output>();
	private LinkedList<Parameter> parameters = new LinkedList<Parameter>();

	public static class Name {
				
		private String id = null;
		private String displayName = null;
		private String prefix;
		private String postfix;

		public static Name createEmptyName() {
			return new Name(null, null, null, null);			
		}

		public static Name createName(String name, String humanReadableName) {
			return new Name(name, null, null, humanReadableName);
		}

		public static Name createName(String name) {
			return new Name(name, null, null, name);
		}

		public static Name createNameSet(String prefix, String postfix, String humanReadableName) {
			return new Name(null, prefix, postfix, humanReadableName);
		}

		public static Name createNameSet(String prefix, String postfix) {
			return new Name(null, prefix, postfix, prefix + SADLSyntax.INPUT_SET_DESIGNATOR + postfix);
		}

		private Name(String id, String prefix, String postfix, String displayName) {
			this.id = id;
			this.prefix = prefix;
			this.postfix = postfix;
			this.displayName = displayName;
		}

		public void setPrefix(String prefix) {
			this.prefix = prefix;
		}

		public void setPostfix(String postfix) {
			this.postfix = postfix;
		}

		public String getID() {
			return id;
		}

		public String getDisplayName() {
			return displayName;
		}

		public void setID(String id) {
			this.id = id;
		}

		public void setDisplayName(String displayName) {
			this.displayName = displayName;
		}
		
		/**
		 * The prefix part of a spliced id of an input set. For a regular input 
		 * returns null.
		 */
		public String getPrefix() {
			return prefix;
		}

		/**
		 * The postfix part of a spliced id of an input set. For a regular input 
		 * returns null.
		 */
		public String getPostfix() {
			return postfix;
		}
		
		/**
		 * Return true iff this is an input set (not a regular input).
		 */
		public boolean isNameSet( ) {
			return id == null;
		}
		
		@Override
		public String toString() {
			
			String firstPart;
			if (isNameSet()) {
				firstPart = getPrefix() + "{...}" + getPostfix();
			} else {
				firstPart = id;
			}
			
			return "\"" + firstPart + "\": \"" + displayName + "\"";
		}


	}
	

	public static class File {
		
		private Name name;
		private boolean optional = false;

		public File(Name name) {
			this.name = name;
		}

		public boolean isOptional() {
			return optional;
		}

		public void setName(Name name) {
			this.name = name;
		}

		/**
		 * The id of a regular input. For an input set returns null.
		 */
		public Name getName() {
			return name;
		}
		
		public void setOptional(boolean optional) {
			this.optional  = optional;
		}		
		
	}

	/**
	 * Input file description.
	 */
	public static class Input extends File {
		
		private InputType type;

		public Input(InputType type, Name name) {
			super(name);
			this.type = type;
		}

		public Input() {
			super(Name.createEmptyName());
		}

		public void setType(InputType type) {
			this.type = type;
		}

		public InputType getType() {
			return type;
		}
	}

	/**
	 * Output file description.
	 */
	public static class Output extends File {
		
		public Output(Name name) {
			super(name);
		}

		public Output() {
			super(Name.createEmptyName());
		}
	}

	/**
	 * Analysis tool parameter description.
	 *
	 */
	public static class Parameter {

		private Name name; 
		private ParameterType type; 
		private String[] selectionOptions; 
		private String from;
		private String to;
		private String defaultValue; 
		private String comment;
		

		public Parameter(Name name, ParameterType type, String[] selectionOptions,
				String from, String to, String defaultValue, String comment) {
			this.name = name;
			this.type = type;
			this.selectionOptions = selectionOptions;
			this.from = from;
			this.to = to;
			this.defaultValue = defaultValue;
			this.comment = comment;
		}
		
		public Name getName() {
			return name;
		}
		
		public ParameterType getType() {
			return type;
		}
		
		public String[] getSelectionOptions() {
			return selectionOptions;
		}
		
		public String getFrom() {
			return from;
		}
		
		public String getTo() {
			return to;
		}
		
		public String getDefaultValue() {
			return defaultValue;
		}
		
		public String getComment() {
			return comment;
		}		
	}

	/**
	 * Returns a new (mostly empty) object presentation for parsed VVSADL. 
	 */
	public SADLDescription(Name name, String comment) {
		super();
		this.name = name;
		this.comment = comment;
	}

	/**
	 * Returns a new (mostly empty) object presentation for parsed VVSADL. 
	 */
	public SADLDescription(Name name, String packageName, String comment) {
		super();
		this.name = name;
		this.comment = comment;
		// skip packageName
	}

	public void addInput(Input input) {
		inputs.add(input);
	}

	public void addMetaInput(Input input) {
		metaInputs.add(input);
	}

	public void addOutput(Output metaOutput) {
		outputs.add(metaOutput);
	}
	
	public void addMetaOutput(Output metaOutput) {
		metaOutputs.add(metaOutput);
	}
	
	public void addParameter(Parameter parameter) {
		parameters.add(parameter);
	}

	public Name getName() {
		return name;
	}

	public String getComment() {
		return comment;
	}

	public List<Input> inputs() {
		return inputs;
	}
	
	public List<Input> metaInputs() {
		return metaInputs;
	}
	
	public List<Output> outputs() {
		return outputs;		
	}
	
	public List<Output> metaOutputs() {
		return metaOutputs;	
	}

	public List<Parameter> parameters() {
		return parameters;
	}

	public void addInputs(List<Input> inputCollection) {
		inputs.addAll(inputCollection);		
	}

	public void addMetaInputs(List<Input> inputCollection) {
		metaInputs.addAll(inputCollection);		
	}

	public void addOutputs(List<Output> outputCollection) {
		outputs.addAll(outputCollection);		
	}
	
	public void addMetaOutputs(List<Output> outputCollection) {
		metaOutputs.addAll(outputCollection);		
	}

	public String getPackageName() {
		// TODO Auto-generated method stub
		return null;
	}
	
	/**
	 * @see SADLGenerator#generate(SADLDescription)
	 */
	@Override
	public String toString() {
		return SADLGenerator.generate(this);
	}
	
}
