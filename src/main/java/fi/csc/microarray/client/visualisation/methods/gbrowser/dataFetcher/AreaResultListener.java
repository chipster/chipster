package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;

/**
 * Listens to results of an area request. Area request are sent from the processing layer 
 * to the view layer.
 * 
 * @author Petri Klemelä
 *
 */
public interface AreaResultListener {

	public void processAreaResult(AreaResult areaResult);
	
}
