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
	private LinkedList<Parameter> parameters = new LinkedList<Parameter>();

	/**
	 * Name for some description entity. Name consists of two major parts:
	 * identifier and display name. The latter should be used when parameter
	 * name is shown to an end-user.
	 * <p>
	 * Name can also contain a prefix and a postfix. This way a single name
	 * can describe a set of names that begin with a certain prefix and at
	 * the same time end with a certain postfix.
	 * 
	 * @author Aleksi Kallio
	 *
	 */
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
		 * Determine if this name defines a set of similar names by using
		 * prefixes and postfixes.
		 * 
		 * @return true if this name defines a set of names, false otherwise.
		 */
		public boolean isSpliced() {
		    return (getPrefix() != null) || (getPostfix() != null); 
		}
		
		/**
		 * Return true iff this is an input set (not a regular input).
		 */
		public boolean isNameSet() {
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
		private boolean isOptional;
		protected String comment;

		public Entity(Name name, boolean isOptional) {
            this.name = name;
            this.isOptional = isOptional;
		}

        public boolean isOptional() {
			return isOptional;
		}

		public void setOptional(boolean isOptional) {
			this.isOptional  = isOptional;
		}		

		public void setName(Name name) {
			this.name = name;
		}

		public Name getName() {
			return name;
		}

		public void setComment(String comment) {
			this.comment = comment;
		}

		public String getComment() {
			return comment;
		}
	}
	
	public static class IOEntity extends Entity {

		private boolean isMeta;
		
		public IOEntity(Name name, boolean isOptional, boolean isMeta) {
			super(name, isOptional);
            this.isMeta = isMeta;
		}

        public boolean isMeta() {
        	return isMeta;
        }

		public void setMeta(boolean isMeta) {
			this.isMeta  = isMeta;
		}		

	}

	/**
	 * Input file description.
	 */
	public static class Input extends IOEntity {
		
		private InputType type;

		public Input() {
			this(null, Name.createEmptyName(), false);
		}

		public Input(InputType type, Name name) {
		    this(type, name, false);
		}

		public Input(InputType type, Name name, boolean optional) {
		    this(type, name, optional, false);
		}


		public Input(InputType type, Name name, boolean optional, boolean isMeta) {
            super(name, optional, isMeta);
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
	public static class Output extends IOEntity {
		
		public Output() {
			this(Name.createEmptyName());
		}
		
		public Output(Name name) {
			this(name, false);
		}

        public Output(Name name, boolean optional) {
            this(name, optional, false);
        }

        public Output(Name name, boolean optional, boolean isMeta) {
            super(name, optional, isMeta);
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

		public Parameter(Name name, ParameterType type, Name[] selectionOptions,
				String from, String to, String defaultValue) {
			this(name, type, selectionOptions, from, to, defaultValue, null);
		}

		public Parameter(Name name, ParameterType type, Name[] selectionOptions,
				String from, String to, String defaultValue, String comment) {
			this(name, type, selectionOptions, from, to,
			     defaultValue == null ? new String[] {} : new String[] {defaultValue},
			     comment);
		}

		public Parameter(Name name, ParameterType type, Name[] selectionOptions,
				String from, String to, String[] defaultValues) {
			this(name, type, selectionOptions, from, to, defaultValues, null);
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
	}

	public SADLDescription(Name name) {
		this(name, null);
	}

	/**
	 * Returns a new (mostly empty) object presentation for parsed SADL. 
	 */
	public SADLDescription(Name name, String comment) {
		super();
		this.name = name;
		this.comment = comment;
	}

	public void addInput(Input input) {
		inputs.add(input);
	}

	public void addOutput(Output metaOutput) {
		outputs.add(metaOutput);
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

	public void setComment(String comment) {
		this.comment = comment;
	}

	public List<Input> inputs() {
		return inputs;
	}
	
	public List<Output> outputs() {
		return outputs;		
	}
	
	public List<Parameter> parameters() {
		return parameters;
	}

	public void addInputs(List<Input> inputCollection) {
		inputs.addAll(inputCollection);		
	}

	public void addOutputs(List<Output> outputCollection) {
		outputs.addAll(outputCollection);		
	}

	/**
	 * @see SADLGenerator#generate(SADLDescription)
	 */
	@Override
	public String toString() {
		return SADLGenerator.generate(this);
	}

	/**
	 * Deprecated
	 */
    // FIXME remove after not using VVSADL anymore
	public void setID(String id) {
		this.name.setID(id);
	}
	
}
