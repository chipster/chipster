package fi.csc.microarray.module.chipster;

import java.util.List;

import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLParser;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLSyntax.ParameterType;

public class ChipsterVVSADLParser extends SADLParser {


	public ChipsterVVSADLParser() {
		this(null);
	}
	
	public ChipsterVVSADLParser(String filename) {
		super(filename);
		addInputType(ChipsterInputTypes.AFFY);
		addInputType(ChipsterInputTypes.CDNA);
		addInputType(ChipsterInputTypes.GENE_EXPRS);
		addInputType(ChipsterInputTypes.GENELIST);
		addInputType(ChipsterInputTypes.PHENODATA);
	}
	
	public static class Validator {
		
		public void validate(String filename, String vvsadl) throws ParseException {
			ChipsterVVSADLParser parser = new ChipsterVVSADLParser(filename);
			List<SADLDescription> descriptions = parser.parseMultiple(vvsadl);
			for (SADLDescription description : descriptions) {
				checkParsedContent(description);
			}
		}
		
		private void checkParsedContent(SADLDescription description) {
						
			for (Parameter parameter : description.parameters()) {
				if (parameter.getType() == ParameterType.ENUM) {
					// check that enum is not empty
					if (parameter.getSelectionOptions() == null || parameter.getSelectionOptions().length == 0) {
						throw new RuntimeException("enum parameter " + parameter.getName() + " has no options");
					}
					// check that enum default value is legal
					if (parameter.getDefaultValue() != null) {
						boolean found = false;
						for (String value : parameter.getSelectionOptions()) {
							if (parameter.getDefaultValue().equals(value)) {
								found = true;
								break;
							}
						}
						if (!found) {
							throw new RuntimeException("enum parameter " + parameter.getName() + " has undefined default value \"" + parameter.getDefaultValue() + "\"");
						}
					}
				}
			}
			
		}
	};
}
