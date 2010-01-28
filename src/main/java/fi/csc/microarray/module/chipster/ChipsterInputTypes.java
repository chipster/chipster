package fi.csc.microarray.module.chipster;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.description.VVSADLSyntax.InputType;
import fi.csc.microarray.exception.MicroarrayException;

public class ChipsterInputTypes {

	
	public static final InputType CDNA = new InputType() {

		public String getName() {			
			return "CDNA";
		}

		public boolean isTypeOf(DataBean dataBean) {
			return dataBean.queryFeatures("/column/sample").exists();
		}

		public boolean isMetadata() {
			return false;
		}		

	};
	
	
	public static final InputType AFFY = new InputType() {

		public String getName() {
			return "AFFY";
		}

		public boolean isTypeOf(DataBean dataBean) {
			return dataBean.isContentTypeCompatitible("application/cel");
		}

		public boolean isMetadata() {
			return false;
		}

	};
	
	public static final InputType GENE_EXPRS = new InputType() {

		public String getName() {
			return "GENE_EXPRS";
		}

		public boolean isTypeOf(DataBean dataBean) {
			try {
				Table chips = dataBean.queryFeatures("/column/chip.*").asTable();
				return chips != null && chips.getColumnCount() > 0;
			} catch (MicroarrayException e) {
				throw new RuntimeException(e);
			}
		}

		public boolean isMetadata() {
			return false;
		}
		
	};

	public static final InputType GENELIST = new InputType() {

		public String getName() {
			return "GENELIST";
		}

		public boolean isTypeOf(DataBean dataBean) {
			return dataBean.queryFeatures("/identifier").exists();
		}

		public boolean isMetadata() {
			return false;
		}
	
	};

	public static final InputType PHENODATA = new InputType() {

		public String getName() {
			return "PHENODATA";
		}

		public boolean isTypeOf(DataBean dataBean) {
			return dataBean.queryFeatures("/phenodata").exists();
		}
		
		public boolean isMetadata() {
			return true;
		}
	};

	public static boolean hasRawType(DataBean data) {
		return AFFY.isTypeOf(data) || CDNA.isTypeOf(data);
	}
}
