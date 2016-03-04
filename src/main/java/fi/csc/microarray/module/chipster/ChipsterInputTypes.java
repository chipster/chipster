package fi.csc.microarray.module.chipster;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.description.SADLSyntax.InputType;
import fi.csc.microarray.exception.MicroarrayException;

public class ChipsterInputTypes {
	
	public static final InputType CDNA = new InputType("CDNA");	
	public static final InputType AFFY = new InputType("AFFY");	
	public static final InputType GENE_EXPRS = new InputType("GENE_EXPRS");	
	public static final InputType GENELIST = new InputType("GENELIST");	
	public static final InputType PHENODATA = new InputType("PHENODATA");	
	public static final InputType BAM = new InputType("BAM");	
	public static final InputType FASTA = new InputType("FASTA");	
	public static final InputType GTF = new InputType("GTF");	
	public static final InputType MOTHUR_OLIGOS = new InputType("MOTHUR_OLIGOS");	
	public static final InputType MOTHUR_NAMES = new InputType("MOTHUR_NAMES");
	public static final InputType MOTHUR_GROUPS = new InputType("MOTHUR_GROUPS");
	
	public static boolean isTypeOf(DataBean dataBean, InputType type) {
		
		if (CDNA.equals(type)) {
			return dataBean.queryFeatures("/column/sample").exists();
		}
		if (AFFY.equals(type)) {
			return dataBean.isContentTypeCompatitible("application/cel");
		}
		if (GENE_EXPRS.equals(type)) {
			try (Table chips = dataBean.queryFeatures("/column/chip.*").asTable()) {
					return chips != null && chips.getColumnCount() > 0;
			} catch (MicroarrayException e) {
				throw new RuntimeException(e);
			}
		}
		if (GENELIST.equals(type)) {
			return dataBean.queryFeatures("/identifier").exists();
		}
		if (PHENODATA.equals(type)) {
			return dataBean.queryFeatures("/phenodata").exists();
		}
		if (BAM.equals(type)) {
			return dataBean.isContentTypeCompatitible("application/bam");
		}
		if (FASTA.equals(type)) {
			return dataBean.isContentTypeCompatitible("chemical/x-fasta");
		}
		if (GTF.equals(type)) {
			return dataBean.isContentTypeCompatitible("text/gtf");
		}
		if (MOTHUR_OLIGOS.equals(type)) {
			return dataBean.isContentTypeCompatitible("text/mothur-oligos");
		}
		if (MOTHUR_NAMES.equals(type)) {
			return dataBean.isContentTypeCompatitible("text/mothur-names");
		}
		if (MOTHUR_GROUPS.equals(type)) {
			return dataBean.isContentTypeCompatitible("text/mothur-groups");
		}
		return false;
}

	
	public static boolean hasRawType(DataBean data) {
		return isTypeOf(data, AFFY) || isTypeOf(data, CDNA);
	}
}
