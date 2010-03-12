package fi.csc.microarray.description;

import java.util.LinkedList;
import java.util.List;

import sun.net.www.http.Hurryable;

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

	private AnnotatedName name;
	private String comment;
	
	private LinkedList<Input> inputs = new LinkedList<Input>();
	private LinkedList<String> outputs = new LinkedList<String>();
	private LinkedList<Input> metaInputs = new LinkedList<Input>();
	private LinkedList<String> metaOutputs = new LinkedList<String>();
	private LinkedList<Parameter> parameters = new LinkedList<Parameter>();

	public static class AnnotatedName {
		private String name = null;
		private String humanReadableName = null;
		private String prefix;
		private String postfix;

		public static AnnotatedName createEmptyName() {
			return new AnnotatedName(null, null, null, null);			
		}

		public static AnnotatedName createName(String name, String humanReadableName) {
			return new AnnotatedName(name, null, null, humanReadableName);
		}
		
		public static AnnotatedName createNameSet(String prefix, String postfix, String humanReadableName) {
			return new AnnotatedName(null, prefix, postfix, humanReadableName);
		}
		
		private AnnotatedName(String name, String prefix, String postfix, String humanReadableName) {
			this.name = name;
			this.prefix = prefix;
			this.postfix = postfix;
			this.humanReadableName = humanReadableName;
		}

		public void setPrefix(String prefix) {
			this.prefix = prefix;
		}

		public void setPostfix(String postfix) {
			this.postfix = postfix;
		}

		public String getName() {
			return name;
		}

		public String getHumanReadableName() {
			return humanReadableName;
		}

		public void setName(String name) {
			this.name = name;
		}

		public void setHumanReadableName(String humanReadableName) {
			this.humanReadableName = humanReadableName;
		}
		
		/**
		 * The prefix part of a spliced name of an input set. For a regular input 
		 * returns null.
		 */
		public String getPrefix() {
			return prefix;
		}

		/**
		 * The postfix part of a spliced name of an input set. For a regular input 
		 * returns null.
		 */
		public String getPostfix() {
			return postfix;
		}
		
		/**
		 * Return true iff this is an input set (not a regular input).
		 */
		public boolean isNameSet( ) {
			return name == null;
		}
		
		@Override
		public String toString() {
			
			String firstPart;
			if (isNameSet()) {
				firstPart = getPrefix() + "{...}" + getPostfix();
			} else {
				firstPart = name;
			}
			
			return "\"" + firstPart + "\": \"" + humanReadableName + "\"";
		}


	}
	
	/**
	 * Input file description.
	 */
	public static class Input {
		
		private InputType type;
		private AnnotatedName name;
		private boolean optional = false;

		public Input(InputType type, AnnotatedName name) {
			this.type = type;
			this.name = name;
		}

		public Input() {
			this.name = AnnotatedName.createEmptyName();
		}

		public boolean isOptional() {
			return optional;
		}

		public void setType(InputType type) {
			this.type = type;
		}

		public void setName(AnnotatedName name) {
			this.name = name;
		}


		/**
		 * The type of this input. 
		 */
		public InputType getType() {
			return type;
		}
		
		/**
		 * The name of a regular input. For an input set returns null.
		 */
		public AnnotatedName getAnnotatedName() {
			return name;
		}
		
		public void setOptional(boolean optional) {
			this.optional  = optional;
		}		
	}
	
	/**
	 * Analysis tool parameter description.
	 *
	 */
	public static class Parameter {

		private AnnotatedName name; 
		private ParameterType type; 
		private String[] selectionOptions; 
		private String from;
		private String to;
		private String defaultValue; 
		private String comment;
		

		public Parameter(AnnotatedName name, ParameterType type, String[] selectionOptions,
				String from, String to, String defaultValue, String comment) {
			this.name = name;
			this.type = type;
			this.selectionOptions = selectionOptions;
			this.from = from;
			this.to = to;
			this.defaultValue = defaultValue;
			this.comment = comment;
		}
		
		public AnnotatedName getAnnotatedName() {
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
	public SADLDescription(AnnotatedName name, String comment) {
		super();
		this.name = name;
		this.comment = comment;
	}

	/**
	 * Returns a new (mostly empty) object presentation for parsed VVSADL. 
	 */
	public SADLDescription(AnnotatedName name, String packageName, String comment) {
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

	public void addOutput(String metaOutput) {
		outputs.add(metaOutput);
	}
	
	public void addMetaOutput(String metaOutput) {
		metaOutputs.add(metaOutput);
	}
	
	public void addParameter(Parameter parameter) {
		parameters.add(parameter);
	}

	public AnnotatedName getAnnotatedName() {
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
	
	public List<String> outputs() {
		return outputs;		
	}
	
	public List<String> metaOutputs() {
		return metaOutputs;	
	}

	public List<Parameter> parameters() {
		return parameters;
	}

	/**
	 * Creates a VVSADL source code representation of this parsed syntax object.
	 * Due to whitespace etc. the returned code might not be identical to the original
	 * source. However if the returned String is used to create a new parsed syntax, it 
	 * should return the exactly same string.
	 * 
	 * @return VVSADL source representation
	 */
	@Override
	public String toString() {
		
		String string =	"TOOL " + this.getAnnotatedName() + " (" + this.getComment() + ")\n";
		
		string += toStringForInputs("INPUT", inputs());		
		string += toStringForInputs("METAINPUT", metaInputs());
		
		string += toStringForOutputs("OUTPUT", outputs());		
		string += toStringForOutputs("METAOUTPUT", metaOutputs());

		if (!parameters().isEmpty()) {
			for (Parameter parameter: parameters()) {
				String paramString = "PARAMETER " + parameter.getAnnotatedName() + " TYPE ";
				
				if (parameter.getType() == ParameterType.ENUM) {
					paramString += "[";
					boolean first = true;
					for (String option : parameter.getSelectionOptions()) {
						if (!first) {
							paramString += ", ";
						} else {
							first = false;
						}
						paramString += option;
					}
					paramString += "] ";
					
				} else {
					paramString += parameter.getType() + " ";

					if (parameter.getFrom() != null) {
						paramString += "FROM " + parameter.getFrom() + " "; 
					}

					if (parameter.getTo() != null) {
						paramString += "TO " + parameter.getTo() + " "; 
					} 
				}	
				
				if (parameter.getDefaultValue() != null) {
					paramString += "DEFAULT " + parameter.getDefaultValue() + " "; 
				}
				
				paramString += "(" + parameter.getComment() + ")";
				
				string += paramString + "\n";
			}			
		}		
		
		return string;
	}

	private String toStringForOutputs(String header, List<String> outputList) {
		String string = "";
		if (!outputList.isEmpty()) {
			for (String output : outputList) {
				string += header + " " + output + "\n";
			}
		}
		return string;
	}

	private String toStringForInputs(String header, List<Input> inputList) {
		String string = "";
		if (!inputList.isEmpty()) {
			for (Input input : inputList) {
				string += header + " " + input.getAnnotatedName().toString() + " TYPE " + input.getType().getName() + "\n";
			}
			
		}
		return string;
	}

	public void addInputs(List<Input> inputCollection) {
		inputs.addAll(inputCollection);		
	}

	public void addMetaInputs(List<Input> inputCollection) {
		metaInputs.addAll(inputCollection);		
	}

	public void addOutputs(List<String> outputCollection) {
		outputs.addAll(outputCollection);		
	}
	
	public void addMetaOutputs(List<String> outputCollection) {
		metaOutputs.addAll(outputCollection);		
	}

	public String getPackageName() {
		// TODO Auto-generated method stub
		return null;
	}
}
