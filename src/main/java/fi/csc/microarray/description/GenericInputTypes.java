package fi.csc.microarray.description;

import fi.csc.microarray.databeans.Dataset;
import fi.csc.microarray.description.SADLSyntax.InputType;

public class GenericInputTypes {

	
	public static final InputType GENERIC = new InputType() {

		public String getName() {			
			return "GENERIC";
		}

		public boolean isTypeOf(Dataset dataBean) {
			return true;
		}

		public boolean isMetadata() {
			return false;
		}
	
	};	
}
