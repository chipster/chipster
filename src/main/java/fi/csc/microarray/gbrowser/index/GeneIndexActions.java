package fi.csc.microarray.gbrowser.index;

/**
 * Gene indexing tools. 
 * 
 * @author Petri Klemelä
 */

import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.QueueManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class GeneIndexActions implements AreaResultListener {

	public interface GeneLocationListener {
		public void geneLocation(Region geneRegion);
	}

	private QueueManager queueManager;
	private LineDataSource gtfDataSource;
	private Map<String, GeneLocationListener> listenerMap = new HashMap<String, GeneLocationListener>();

	public GeneIndexActions(QueueManager queueManager, LineDataSource gtfDataSource) {

		this.queueManager = queueManager;
		this.gtfDataSource = gtfDataSource;

		queueManager.addResultListener(gtfDataSource, this);

	}

	/**
	 * getting location of a gene
	 * @param chr 
	 * @throws SQLException 
	 */
	public void requestLocation(String gene, Chromosome chr, GeneLocationListener listener) {

		queueManager.addAreaRequest(gtfDataSource, new GeneRequest(gene, chr), false);

		listenerMap.put(gene, listener);

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

			listenerMap.get(geneResult.getSearchString()).geneLocation(geneResult.getGeneLocation());
			listenerMap.remove(geneResult.getSearchString());
		}
	}
}
