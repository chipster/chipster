package fi.csc.microarray.description;

import java.util.List;

import fi.csc.microarray.description.SADLDescription.Entity;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLDescription.Output;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLSyntax.ParameterType;
import fi.csc.microarray.util.Strings;

/**
 * Generates SADL source code from parsed objects.
 * 
 * @author Aleksi Kallio
 *
 */
public class SADLGenerator {

	/**
	 * Creates a SADL source code representation of parsed syntax object (SADLDescription).
	 * Due to whitespace etc. the returned code might not be identical to the original
	 * source. However if the returned String is used to create a new parsed syntax, it 
	 * should return the exactly same string.
	 * 
	 * @return SADL source representation
	 */
	public static String generate(SADLDescription sadl) {
		
		String string =	"TOOL " + quoteIfNeeded(sadl.getCategory()) + " / " + generateName(sadl.getName()) + " (" + escapeIfNeeded(sadl.getComment()) + ")\n";
		
		string += generateInputs("INPUT", sadl.inputs());		
		string += generateInputs("METAINPUT", sadl.metaInputs());
		
		string += generateOutputs("OUTPUT", sadl.outputs());		
		string += generateOutputs("METAOUTPUT", sadl.metaOutputs());

		if (!sadl.parameters().isEmpty()) {
			for (Parameter parameter: sadl.parameters()) {
				String paramString = "PARAMETER " + generateOptional(parameter) + parameter.getName() + " TYPE ";
				
				if (parameter.getType() == ParameterType.ENUM) {
					paramString += "[";
					boolean first = true;
					for (Name option : parameter.getSelectionOptions()) {
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
				}
				
				if (parameter.getFrom() != null) {
					paramString += "FROM " + parameter.getFrom() + " "; 
				}

				if (parameter.getTo() != null) {
					paramString += "TO " + parameter.getTo() + " "; 
				} 

				if (parameter.getDefaultValues().length > 0) {
					paramString += "DEFAULT ";
					boolean first = true;
					for (String defaultValue : parameter.getDefaultValues()) {
						paramString += first ? "" : ",";
						paramString += quoteIfNeeded(defaultValue) + " ";
						first = false;
					}
				}

				paramString += "(" + escapeIfNeeded(parameter.getComment()) + ")";
				
				string += paramString + "\n";
			}			
		}		
		
		return string;
	}

	private static String generateOutputs(String header, List<Output> outputList) {
		String string = "";
		if (!outputList.isEmpty()) {
			for (Output output : outputList) {
				string += header + " " + generateOptional(output) + generateName(output.getName()) +  "\n";
			}
		}
		return string;
	}

	private static String generateInputs(String header, List<Input> inputList) {
		String string = "";
		if (!inputList.isEmpty()) {
			for (Input input : inputList) {
				string += header + " " + generateOptional(input) + generateName(input.getName()) + " TYPE " + input.getType().getName() + "\n";
			}
			
		}
		return string;
	}
	
	private static String generateOptional(Entity entity) {
		return entity.isOptional() ? "OPTIONAL " : "";		
	}
	
	public static String generateName(Name name) {
		
		String firstPart;
		if (name.isNameSet()) {
			firstPart = name.getPrefix() + "{...}" + name.getPostfix();
		} else {
			firstPart = name.getID();
		}
		
		System.out.println(name.getID());
		return quoteIfNeeded(firstPart) + ": " + quoteIfNeeded(name.getDisplayName());
	}
	
	private static String quoteIfNeeded(String string) {
		System.out.println(string);
		if (string.contains(" ") || Strings.containsAnyOf(string, true, SADLTokeniser.tokenEndingOperators())) {
			return "\"" + escapeIfNeeded(string) + "\"";
		} else {
			return string;
		}
	}
	
	private static String escapeIfNeeded(String string) {
		for (String operator : SADLTokeniser.blockEndingOperators()) {
			string = string.replace(operator, SADLSyntax.ESCAPE + operator);
		}
		
		return string;
	}


}
