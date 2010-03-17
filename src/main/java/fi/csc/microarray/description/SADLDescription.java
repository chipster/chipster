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
			
			return hyphenate(firstPart) + ": " + hyphenate(displayName);
		}
		
		private String hyphenate(String name) {
			if (name.contains(" ")) {
				return "\"" + name + "\"";
			} else {
				return name;
			}
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
			super(name, false);
		}

		public Output() {
			super(Name.createEmptyName(), false);
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
		private String defaultValue; 
		private String comment;
		

		public Parameter(Name name, ParameterType type, Name[] selectionOptions,
				String from, String to, String defaultValue, String comment) {
			super(name, false);
			this.type = type;
			this.selectionOptions = selectionOptions;
			this.from = from;
			this.to = to;
			this.defaultValue = defaultValue;
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
			return defaultValue;
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

}
