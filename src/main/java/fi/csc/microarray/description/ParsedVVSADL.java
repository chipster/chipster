package fi.csc.microarray.description;

import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.description.VVSADLSyntax.InputType;
import fi.csc.microarray.description.VVSADLSyntax.ParameterType;

public class ParsedVVSADL {

	private String name;
	private String packageName; 
	private String comment;
	
	private LinkedList<Input> inputs = new LinkedList<Input>();
	private LinkedList<String> outputs = new LinkedList<String>();
	private LinkedList<Input> metaInputs = new LinkedList<Input>();
	private LinkedList<String> metaOutputs = new LinkedList<String>();
	private LinkedList<Parameter> parameters = new LinkedList<Parameter>();
	
	public static class Input {
		
		private InputType type;
		private String name;
		private String prefix;
		private String postfix;
		
		public static Input createInput(InputType type, String name) {
			return new Input(type, name, null, null);
		}

		public static Input createInputSet(InputType type, String prefix, String postfix) {
			return new Input(type, null, prefix, postfix);
		}

		private Input(InputType type, String name, String prefix, String postfix) {
			this.type = type;
			this.name = name;
			this.prefix = prefix;
			this.postfix = postfix;
		}

		public InputType getType() {
			return type;
		}
		
		public String getName() {
			return name;
		}
		
		public String getPrefix() {
			return prefix;
		}
		public String getPostfix() {
			return postfix;
		}
		public boolean isInputSet( ) {
			return name == null;
		}
	}
	
	public static class Parameter {

		private String name; 
		private ParameterType type; 
		private String[] selectionOptions; 
		private String from;
		private String to;
		private String defaultValue; 
		private String comment;
		

		public Parameter(String name, ParameterType type, String[] selectionOptions,
				String from, String to, String defaultValue, String comment) {
			this.name = name;
			this.type = type;
			this.selectionOptions = selectionOptions;
			this.from = from;
			this.to = to;
			this.defaultValue = defaultValue;
			this.comment = comment;
		}
		
		public String getName() {
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

	public ParsedVVSADL(String name, String packageName, String comment) {
		super();
		this.name = name;
		this.packageName = packageName;
		this.comment = comment;
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

	public String getName() {
		return name;
	}

	public String getPackageName() {
		return packageName;
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
		
		String string =	"ANALYSIS \"" + this.getPackageName() + "\"/\"" + this.getName() + "\" (" + this.getComment() + ")\n";
		
		string += toStringForInputs("INPUT", inputs());		
		string += toStringForInputs("METAINPUT", metaInputs());
		
		string += toStringForOutputs("OUTPUT", outputs());		
		string += toStringForOutputs("METAOUTPUT", metaOutputs());

		if (!parameters().isEmpty()) {
			for (Parameter parameter: parameters()) {
				String paramString = "PARAMETER " + parameter.getName() + " ";
				
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
			String outputString = header + " ";
			boolean first = true;
			for (String output : outputList) {
				if (!first) {
					outputString += ", ";
				} else {
					first = false;
				}
				outputString += output;
			}
			
			string += outputString + "\n";
		}
		return string;
	}

	private String toStringForInputs(String header, List<Input> inputList) {
		String string = "";
		if (!inputList.isEmpty()) {
			String inputString = header + " ";
			boolean first = true;
			for (Input input : inputList) {
				if (!first) {
					inputString += ", ";
				} else {
					first = false;
				}
				inputString += input.getType().getName() + " ";
				if (input.isInputSet()) {
					inputString += input.getPrefix() + "[...]" + input.getPostfix();
				} else {
					inputString += input.getName();
				}
			}
			
			string += inputString + "\n";
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
}
