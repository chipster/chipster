package fi.csc.microarray.client.visualisation.methods.gbrowser;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;

public interface RegionListener {
	
	public void RegionChanged(BpCoordRegion bpRegion);
}
