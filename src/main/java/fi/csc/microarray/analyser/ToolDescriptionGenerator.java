package fi.csc.microarray.analyser;

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
	public ToolDescription generate(SADLDescription source, AnalysisHandler analysisHandler) {
		ToolDescription description = new ToolDescription(analysisHandler);
		
		description.setID(source.getName().getID());
		description.setDisplayName(source.getName().getDisplayName());
		description.setComment(source.getComment());

		// not interested in inputs, they were figured out when job was submitted
		// I'm interested in inputs in java jobs
		for (Input input : source.inputs()) {
			description.addInputFile(input.getName().getID());
		}
		
		for (Output output : source.outputs()) {
			description.addOutputFile(output.getName(), output.isOptional());
		}
		
		for (Parameter parameter : source.parameters()) {
			description.addParameter(new ToolDescription.ParameterDescription(parameter.getName().getID(), parameter.getComment(), parameter.getType().isNumeric()));
		}
		
		return description;
	}
}