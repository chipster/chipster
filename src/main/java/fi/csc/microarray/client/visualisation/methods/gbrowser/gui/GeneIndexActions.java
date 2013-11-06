package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

/**
 * Gene indexing tools. 
 * 
 * @author Petri Klemelä
 */

import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;

/**
 * This class does gene search in two steps: first use custom gene-chr file to find the chromosome 
 * of the searched gene and then find the exact location from the Gtf file. The first step is needed 
 * to keep the amount of Gtf-reading reasonably small.
 * 
 * 
 * @author klemela
 *
 */
public class GeneIndexActions implements DataResultListener {

	public interface GeneLocationListener {
		public void geneLocation(Region geneRegion);
	}

	private QueueManager queueManager;
	private DataThread gtfDataSource;
	private DataThread geneDataSource;
	private Map<String, GeneLocationListener> listenerMap = new HashMap<String, GeneLocationListener>();

	public GeneIndexActions(QueueManager queueManager, DataThread gtfDataSource, DataThread geneDataSource) {

		this.queueManager = queueManager;
		this.gtfDataSource = gtfDataSource;
		this.geneDataSource = geneDataSource;

		initializeDataResultListeners();
	}

	/**
	 * getting location of a gene
	 * @throws SQLException 
	 */
	public void requestLocation(String gene, GeneLocationListener listener) {

		// Listener is called later when we have both two search results
		listenerMap.put(gene, listener);
		
		// Search for the chromosome of the gene
		queueManager.addDataRequest(geneDataSource, new GeneRequest(gene, null), null);
	}


	private void requestLocation(String gene, Chromosome chr) {
		
		// We know the chromosome, but search for location of the gene
		queueManager.addDataRequest(gtfDataSource, new GeneRequest(gene, chr), null);
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
	public void processDataResult(DataResult dataResult) {
		if (dataResult instanceof GeneResult) {
			GeneResult geneResult = (GeneResult) dataResult;
			
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

	public void initializeDataResultListeners() {
		queueManager.addDataResultListener(gtfDataSource, this);
		queueManager.addDataResultListener(geneDataSource, this);
	}
}
