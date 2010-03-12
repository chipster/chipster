package fi.csc.microarray.description;

import java.util.List;

import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLSyntax.ParameterType;

public class SADLGenerator {

	/**
	 * Creates a VVSADL source code representation of description parsed syntax object.
	 * Due to whitespace etc. the returned code might not be identical to the original
	 * source. However if the returned String is used to create a new parsed syntax, it 
	 * should return the exactly same string.
	 * 
	 * @return VVSADL source representation
	 */
	public String generate(SADLDescription description) {
		
		String string =	"ANALYSIS \"" + description.getPackageName() + "\"/\"" + description.getAnnotatedName() + "\" (" + description.getComment() + ")\n";
		
		string += generateInputs("INPUT", description.inputs());		
		string += generateInputs("METAINPUT", description.metaInputs());
		
		string += generateOutputs("OUTPUT", description.outputs());		
		string += generateOutputs("METAOUTPUT", description.metaOutputs());

		if (!description.parameters().isEmpty()) {
			for (Parameter parameter: description.parameters()) {
				String paramString = "PARAMETER " + parameter.getAnnotatedName() + " ";
				
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

	private String generateOutputs(String header, List<String> outputList) {
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

	private String generateInputs(String header, List<Input> inputList) {
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
				if (input.getAnnotatedName().isNameSet()) {
					inputString += input.getAnnotatedName().getPrefix() + "[...]" + input.getAnnotatedName().getPostfix();
				} else {
					inputString += input.getAnnotatedName();
				}
			}
			
			string += inputString + "\n";
		}
		return string;
	}
}
