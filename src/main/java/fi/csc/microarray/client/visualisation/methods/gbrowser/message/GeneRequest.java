package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.HashSet;


/**
 * Class to disguise gene search as a dataRequest. Note that as this really is also DataRequest, so it's
 * easily deleted when the queues are cleared or when requests are combined in QueueManager. This is a problem
 * for user interface, which waits for the search result. Avoid redrawing screen when the gene search is sent to
 * avoid other DataRequest causing these problems.
 * 
 * @author klemela
 *
 */
public class GeneRequest extends DataRequest {
	
	private String gene = null;

	public GeneRequest(String gene, Chromosome chr) {
		super(
				//chr is used only to optimize Gtf reading, every chromosome is search for the gene
				new Region(0l, Long.MAX_VALUE, chr), 
				new HashSet<DataType>(), 
				new DataStatus());
		
		this.gene = gene;
	}
	
	public String getSearchString() {
		return gene;
	}
}
