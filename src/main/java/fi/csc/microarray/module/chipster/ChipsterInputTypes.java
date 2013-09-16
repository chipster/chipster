package fi.csc.microarray.module.chipster;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.description.SADLSyntax.InputType;
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
	
	public static final InputType BAM = new InputType() {

		public String getName() {
			return "BAM";
		}

		public boolean isTypeOf(DataBean dataBean) {
			return dataBean.isContentTypeCompatitible("application/bam");
		}

		public boolean isMetadata() {
			return false;
		}

	};
	
	public static final InputType FASTA = new InputType() {

		public String getName() {
			return "FASTA";
		}

		public boolean isTypeOf(DataBean dataBean) {
			return dataBean.isContentTypeCompatitible("chemical/x-fasta");
		}

		public boolean isMetadata() {
			return false;
		}

	};

	public static final InputType GTF = new InputType() {

		public String getName() {
			return "GTF";
		}

		public boolean isTypeOf(DataBean dataBean) {
			return dataBean.isContentTypeCompatitible("text/gtf");
		}

		public boolean isMetadata() {
			return false;
		}

	};

	public static final InputType MOTHUR_OLIGOS = new InputType() {

		public String getName() {
			return "MOTHUR_OLIGOS";
		}

		public boolean isTypeOf(DataBean dataBean) {
			return dataBean.isContentTypeCompatitible("text/mothur-oligos");
		}

		public boolean isMetadata() {
			return false;
		}

	};

	public static final InputType MOTHUR_NAMES = new InputType() {

		public String getName() {
			return "MOTHUR_NAMES";
		}

		public boolean isTypeOf(DataBean dataBean) {
			return dataBean.isContentTypeCompatitible("text/mothur-names");
		}

		public boolean isMetadata() {
			return false;
		}

	};

	public static final InputType MOTHUR_GROUPS = new InputType() {

		public String getName() {
			return "MOTHUR_GROUPS";
		}

		public boolean isTypeOf(DataBean dataBean) {
			return dataBean.isContentTypeCompatitible("text/mothur-groups");
		}

		public boolean isMetadata() {
			return false;
		}

	};

	
	
	
	public static boolean hasRawType(DataBean data) {
		return AFFY.isTypeOf(data) || CDNA.isTypeOf(data);
	}
}
