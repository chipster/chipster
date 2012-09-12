package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.HashSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;

/**
 * Class to disguise gene search as a areaRequest. Note that as this really is also AreaRequest, so it's
 * easily deleted when the queues are cleared or when requests are combined in QueueManager. This is a problem
 * for user interface, which waits for the search result. Avoid redrawing screen when the gene search is sent to
 * avoid other AreaRequest causing these problems.
 * 
 * @author klemela
 *
 */
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
