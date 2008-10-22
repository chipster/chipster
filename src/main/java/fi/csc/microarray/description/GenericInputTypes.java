package fi.csc.microarray.description;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.description.VVSADLSyntax.InputType;

public class GenericInputTypes {

	
	public static final InputType GENERIC = new InputType() {

		public String getName() {			
			return "GENERIC";
		}

		public boolean isTypeOf(DataBean dataBean) {
			return true;
		}

		public boolean isMetadata() {
			return false;
		}
	
	};	
}
