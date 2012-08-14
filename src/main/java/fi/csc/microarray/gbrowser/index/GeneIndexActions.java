package fi.csc.microarray.gbrowser.index;

/**
 * Gene indexing tools. 
 * 
 * @author Petri Klemelä
 */

import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.QueueManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

/**
 * This class does gene search in two steps: first use custom gene-chr file to find the chromosome 
 * of the searched gene and then find the exact location from the Gtf file. The first step is needed 
 * to keep the amount of Gtf-reading reasonably small.
 * 
 * 
 * @author klemela
 *
 */
public class GeneIndexActions implements AreaResultListener {

	public interface GeneLocationListener {
		public void geneLocation(Region geneRegion);
	}

	private QueueManager queueManager;
	private TabixDataSource gtfDataSource;
	private LineDataSource geneDataSource;
	private Map<String, GeneLocationListener> listenerMap = new HashMap<String, GeneLocationListener>();

	public GeneIndexActions(QueueManager queueManager, TabixDataSource gtfDataSource, LineDataSource geneDataSource) {

		this.queueManager = queueManager;
		this.gtfDataSource = gtfDataSource;
		this.geneDataSource = geneDataSource;

		queueManager.addResultListener(gtfDataSource, this);
		queueManager.addResultListener(geneDataSource, this);
	}

	/**
	 * getting location of a gene
	 * @throws SQLException 
	 */
	public void requestLocation(String gene, GeneLocationListener listener) {

		// Listener is called later when we have both two search results
		listenerMap.put(gene, listener);
		
		// Search for the chromosome of the gene
		queueManager.addAreaRequest(geneDataSource, new GeneRequest(gene, null), false);
	}


	private void requestLocation(String gene, Chromosome chr) {
		
		// We know the chromosome, but search for location of the gene
		queueManager.addAreaRequest(gtfDataSource, new GeneRequest(gene, chr), false);
	}

	public static boolean checkIfNumber(String name) {
		try {
			Integer.parseInt(name);
			return true;
		} catch (NumberFormatException e) {
			try {
				Long.parseLong(name);
				return true;
			} catch (NumberFormatException e1) {
				return false;
			}
		}
	}


	@Override
	public void processAreaResult(AreaResult areaResult) {
		if (areaResult instanceof GeneResult) {
			GeneResult geneResult = (GeneResult) areaResult;
			
			if (geneResult.getGeneLocation() == null) {
				// There isn't such gene, return after first search
				listenerMap.get(geneResult.getSearchString()).geneLocation(geneResult.getGeneLocation());
				listenerMap.remove(geneResult.getSearchString());
			} else if (geneResult.getGeneLocation().start.bp != null) {
				// We have the location, both two searches are done
				listenerMap.get(geneResult.getSearchString()).geneLocation(geneResult.getGeneLocation());
				listenerMap.remove(geneResult.getSearchString());
			} else {
				// We have the chromosome, but no location, search for the latter
				requestLocation(geneResult.getSearchString(), geneResult.getGeneLocation().start.chr);
			}
		}
	}
}
