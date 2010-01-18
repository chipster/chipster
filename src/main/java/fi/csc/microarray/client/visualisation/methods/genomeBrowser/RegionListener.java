package fi.csc.microarray.client.visualisation.methods.genomeBrowser;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.BpCoordRegion;

public interface RegionListener {
	public void RegionChanged(BpCoordRegion bpRegion);
}
