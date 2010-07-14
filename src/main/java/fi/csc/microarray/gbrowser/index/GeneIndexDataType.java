package fi.csc.microarray.gbrowser.index;

//TODO if there are chromosomes with X,Y values, then chromosome type should be changed
public class GeneIndexDataType {
	public Long chromosome = null;
	public Long bpstart = null;
	public Long bpend = null;
	
	public GeneIndexDataType(Long chromosome, Long bpstart, Long bpend){
		this.chromosome = chromosome;
		this.bpstart = bpstart;
		this.bpend = bpend;
	}
	
	public GeneIndexDataType(){
		
	}
}
