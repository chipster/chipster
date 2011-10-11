package fi.csc.microarray.client.visualisation.methods.gbrowser;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;

/**
 * Listens to region change events. Region is changed when Genome browser moves (pan or zoom).
 * 
 * @author Petri Klemelä
 *
 */
public interface RegionListener {
	
	public void regionChanged(BpCoordRegion bpRegion);
}
