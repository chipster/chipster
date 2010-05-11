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
	private String category;

	public static class Name {
				
		private String id = null;
		private String displayName = null;
		private String prefix;
		private String postfix;

		public static Name createEmptyName() {
			return new Name(null, null, null, null);			
		}

		public static Name createName(String name, String displayName) {
			return new Name(name, null, null, displayName);
		}

		public static Name createName(String name) {
			return new Name(name, null, null, name);
		}

		public static Name createNameSet(String prefix, String postfix, String displayName) {
			return new Name(null, prefix, postfix, displayName);
		}

		public static Name createNameSet(String prefix, String postfix) {
			return new Name(null, prefix, postfix, prefix + SADLSyntax.NAME_SET_DESIGNATOR + postfix);
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
		
		/**
		 * @see SADLGenerator#generateName(Name)
		 */
		@Override
		public String toString() {
			return SADLGenerator.generateName(this);
		}
		
	}
	

	public static class Entity {
		
		private Name name;
		private boolean optional;

		public Entity(Name name, boolean optional) {
            this.optional = optional;
            this.name = name;
        }

        public boolean isOptional() {
			return optional;
		}

		public void setName(Name name) {
			this.name = name;
		}

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
	public static class Input extends Entity {
		
		private InputType type;

		public Input(InputType type, Name name) {
		    this(type, name, false);
		}

		public Input() {
			super(Name.createEmptyName(), false);
		}

		public Input(InputType type, Name name, boolean optional) {
            super(name, optional);
            this.type = type;
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
	public static class Output extends Entity {
		
		public Output(Name name) {
			this(name, false);
		}

		public Output() {
			this(Name.createEmptyName(), false);
		}
		
        public Output(Name name, boolean optional) {
            super(name, optional);
        }

	}

	/**
	 * Analysis tool parameter description.
	 *
	 */
	public static class Parameter extends Entity {

		private ParameterType type; 
		private Name[] selectionOptions; 
		private String from;
		private String to;
		private String[] defaultValues; 
		private String comment;
		

		public Parameter(Name name, ParameterType type, Name[] selectionOptions,
				String from, String to, String defaultValue, String comment) {
			this(name, type, selectionOptions, from, to,
			     defaultValue == null ? new String[] {} : new String[] {defaultValue},
			     comment);
		}

		public Parameter(Name name, ParameterType type, Name[] selectionOptions,
				String from, String to, String[] defaultValues, String comment) {
			super(name, false);
			this.type = type;
			this.selectionOptions = selectionOptions;
			this.from = from;
			this.to = to;
			this.defaultValues = defaultValues;
			this.comment = comment;
		}
		
		public ParameterType getType() {
			return type;
		}
		
		public Name[] getSelectionOptions() {
			return selectionOptions;
		}
		
		public String getFrom() {
			return from;
		}
		
		public String getTo() {
			return to;
		}
		
		public String getDefaultValue() {
			if (defaultValues.length != 1) {
				throw new IllegalStateException("there needs to be 1 default value, not " + defaultValues.length);
			}
			return defaultValues[0];
		}

		public String[] getDefaultValues() {
			return defaultValues;
		}

		public String getComment() {
			return comment;
		}		
	}

	/**
	 * Returns a new (mostly empty) object presentation for parsed SADL. 
	 */
	public SADLDescription(Name name, String category, String comment) {
		super();
		this.name = name;
		this.comment = comment;
		this.category = category;
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

	public String getCategory() {
		return this.category;
	}
	
	/**
	 * @see SADLGenerator#generate(SADLDescription)
	 */
	@Override
	public String toString() {
		return SADLGenerator.generate(this);
	}

	public String toStringVerbose() {
		String s = "";
		s += "-------------- sadl description --------------\n";
		s += this.getName().getID() + "\n";
		s += this.getName().getDisplayName() + "\n";
		s += this.getCategory() + "\n";
		
		for (SADLDescription.Input input: this.inputs()) {
			String inputID = input.getName().getID();
			if (inputID == null) {
				inputID = input.getName().getPrefix() + input.getName().getPostfix();
			}
			
			s += inputID + ", " + input.getName().getDisplayName() + ", " + input.getType().getName() + "\n";
			
		}
		
		for (SADLDescription.Output output: this.outputs()) {
			String outputID = output.getName().getID();
			if (outputID == null) {
				outputID = output.getName().getPrefix() + output.getName().getPostfix();
			}

			s += outputID + ", " + output.getName().getDisplayName() + ", " + output.isOptional() + "\n";
		}
		for (SADLDescription.Parameter parameter: this.parameters()) {
			s += parameter.getName().getID() + ", " + parameter.getName().getDisplayName() + ", " + parameter.getType() + "\n";
		}
		
		s += "-------------- sadl description --------------\n";
		return s;
	}
	
}
