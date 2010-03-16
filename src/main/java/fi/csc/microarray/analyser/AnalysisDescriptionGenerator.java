package fi.csc.microarray.analyser;

import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Output;
import fi.csc.microarray.description.SADLDescription.Parameter;

/**
 *  
 * Generator for AnalysisDescription objects. AnalysisDescription objects are
 * compute service specific versions of analysis tools descriptions.
 * 
 * @author Aleksi Kallio
 *
 */
public class AnalysisDescriptionGenerator {

	/**
	 * Converts generic SADLDescription to AnalysisDescription.
	 */
	public AnalysisDescription generate(SADLDescription source, AnalysisHandler analysisHandler) {
		AnalysisDescription description = new AnalysisDescription(analysisHandler);
		
		description.setName(source.getName().getID());
		description.setComment(source.getComment());
		description.setCategory(source.getCategory());

		// not interested in inputs, they were figured out when job was submitted

		for (Output output : source.outputs()) {
			description.addOutputFile(output.getName().getID());
		}
		
		for (Parameter parameter : source.parameters()) {
			description.addParameter(new AnalysisDescription.ParameterDescription(parameter.getName().getID(), parameter.getComment(), parameter.getType().isNumeric()));
		}
		
		return description;
	}

}