package fi.csc.microarray.analyser;

import fi.csc.microarray.description.ParsedVVSADL;
import fi.csc.microarray.description.ParsedVVSADL.Parameter;

/**
 * Listens to parse process and constructs description.
 */
public class AnalysisDescriptionGenerator {

	public AnalysisDescription generate(ParsedVVSADL source, AnalysisHandler analysisHandler) {
		AnalysisDescription description = new AnalysisDescription(analysisHandler);
		
		description.setName(source.getName());
		description.setComment(source.getComment());
		description.setCategory(source.getPackageName());

		// not interested in inputs, they were figured out when job was submitted

		for (String output : source.outputs()) {
			description.addOutputFile(output);
		}
		
		for (Parameter parameter : source.parameters()) {
			description.addParameter(new AnalysisDescription.ParameterDescription(parameter.getName(), parameter.getComment(), parameter.getType().isNumeric()));
		}
		
		return description;
	}

}