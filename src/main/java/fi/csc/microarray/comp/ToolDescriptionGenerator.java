package fi.csc.microarray.comp;

import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Output;
import fi.csc.microarray.description.SADLDescription.Parameter;

/**
 *  
 * Generator for ToolDescription objects. ToolDescription objects are
 * compute service specific versions of analysis tools descriptions.
 * 
 * @author Aleksi Kallio
 *
 */
public class ToolDescriptionGenerator {

	/**
	 * Converts generic SADLDescription to ToolDescription.
	 */
	public ToolDescription generate(SADLDescription source) {
		ToolDescription description = new ToolDescription();
		
		description.setID(source.getName().getID());
		description.setDisplayName(source.getName().getDisplayName());
		description.setComment(source.getDescription());

		// not interested in inputs, they were figured out when job was submitted
		// I'm interested in inputs in java jobs
		for (Input input : source.getInputs()) {
			description.addInputFile(input.getName().getID());
		}
		
		for (Output output : source.getOutputs()) {
			description.addOutputFile(output.getName(), output.isOptional());
		}
		
		for (Parameter parameter : source.getParameters()) {
			description.addParameter(new ToolDescription.ParameterDescription(parameter.getName().getID(), parameter.getDescription(), parameter.getType().isNumeric()));
		}
		
		return description;
	}
}