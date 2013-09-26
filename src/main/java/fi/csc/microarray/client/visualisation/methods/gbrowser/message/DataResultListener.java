package fi.csc.microarray.client.visualisation.methods.gbrowser.message;


/**
 * Listens to results of an data request. Data requests are sent from the processing layer 
 * to the view layer.
 * 
 * @author Petri Klemelä
 *
 */
public interface DataResultListener {

	public void processDataResult(DataResult dataResult);
	
}
