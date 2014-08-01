package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;

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
		
		if (gtfDataSource == null && geneDataSource == null) {
			throw new IllegalArgumentException("Both gene search data sources cannot be null");
		}

		initializeDataResultListeners();
	}

	/**
	 * getting location of a gene
	 * @throws SQLException 
	 */
	public void requestLocation(String gene, GeneLocationListener listener) {

		// Listener is called later when we have both two search results
		listenerMap.put(gene, listener);
		
		if (geneDataSource != null) {
			queueManager.addDataRequest(geneDataSource, new GeneRequest(gene, null), null);
		} else if (gtfDataSource != null) {
			queueManager.addDataRequest(gtfDataSource, new GeneRequest(gene, null), null);
		}
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
			
			listenerMap.get(geneResult.getSearchString()).geneLocation(geneResult.getGeneLocation());
			listenerMap.remove(geneResult.getSearchString());
		}
	}

	public void initializeDataResultListeners() {
		if (gtfDataSource != null) {
			queueManager.addDataResultListener(gtfDataSource, this);
		}
		if (geneDataSource != null) {
			queueManager.addDataResultListener(geneDataSource, this);
		}
	}
}
