package fi.csc.microarray.module.chipster;

import java.util.List;

import fi.csc.microarray.description.ParsedVVSADL;
import fi.csc.microarray.description.VVSADLParser;
import fi.csc.microarray.description.ParsedVVSADL.Parameter;
import fi.csc.microarray.description.VVSADLSyntax.ParameterType;

public class ChipsterVVSADLParser extends VVSADLParser {


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
			List<ParsedVVSADL> descriptions = parser.parseMultiple(vvsadl);
			for (ParsedVVSADL description : descriptions) {
				checkParsedContent(description);
			}
		}
		
		private void checkParsedContent(ParsedVVSADL description) {
						
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
