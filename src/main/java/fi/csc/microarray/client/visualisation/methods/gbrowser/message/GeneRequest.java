package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.HashSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;

public class GeneRequest extends AreaRequest {
	
	private String gene = null;

	public GeneRequest(String gene, Chromosome chr) {
		super(
				//chr is used only to optimize Gtf reading, every chromosome is search for the gene
				new Region(0l, Long.MAX_VALUE, chr), 
				new HashSet<ColumnType>(), 
				new FsfStatus());
		
		this.gene = gene;
	}
	
	public String getSearchString() {
		return gene;
	}
}
