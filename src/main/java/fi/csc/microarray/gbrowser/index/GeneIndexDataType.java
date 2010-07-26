package fi.csc.microarray.gbrowser.index;

/**
 * This is for GeneIndexActions, for saving gene's location 
 */

public class GeneIndexDataType {
	public String chromosome = null;
	public Long bpstart = null;
	public Long bpend = null;
	
	public GeneIndexDataType(String chromosome, Long bpstart, Long bpend){
		this.chromosome = chromosome;
		this.bpstart = bpstart;
		this.bpend = bpend;
	}
	
	public GeneIndexDataType(){
		
	}
}
