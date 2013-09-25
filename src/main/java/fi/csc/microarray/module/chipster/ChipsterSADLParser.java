package fi.csc.microarray.module.chipster;

import java.util.List;

import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLParser;
import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLSyntax.ParameterType;

public class ChipsterSADLParser extends SADLParser {


	public ChipsterSADLParser() {
		this(null);
	}
	
	public ChipsterSADLParser(String filename) {
		super(filename);
		addInputType(ChipsterInputTypes.AFFY);
		addInputType(ChipsterInputTypes.CDNA);
		addInputType(ChipsterInputTypes.GENE_EXPRS);
		addInputType(ChipsterInputTypes.GENELIST);
		addInputType(ChipsterInputTypes.PHENODATA);
		addInputType(ChipsterInputTypes.BAM);
		addInputType(ChipsterInputTypes.FASTA);
		addInputType(ChipsterInputTypes.GTF);
		addInputType(ChipsterInputTypes.MOTHUR_OLIGOS);
		addInputType(ChipsterInputTypes.MOTHUR_NAMES);
		addInputType(ChipsterInputTypes.MOTHUR_GROUPS);
	}
	
	public static class Validator {
		
		public void validate(String filename, String sadl) throws ParseException {
			ChipsterSADLParser parser = new ChipsterSADLParser(filename);
			List<SADLDescription> descriptions = parser.parseMultiple(sadl);
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
					for (String defaultValue : parameter.getDefaultValues()) {
						boolean found = false;
						for (Name value : parameter.getSelectionOptions()) {
							if (defaultValue.equals(value.getID())) {
								found = true;
								break;
							}
						}
						if (!found) {
							throw new RuntimeException("enum parameter " + parameter.getName() + " has undefined default value \"" + defaultValue + "\"");
						}
					}
					
				} else {
					// check that non-enum values do not have multiple default values
					if (parameter.getDefaultValues().length > 1) {
						throw new RuntimeException("non-enum parameter " + parameter.getName() + " has multiple default values");
					}
				}
			}
			
		}
	};
}
