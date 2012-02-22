package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.FileNotFoundException;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.CytobandHandlerThread;

public class CytobandDataSource extends LineDataSource {
	
	private LineDataSource regionDataSource;

	public CytobandDataSource(URL cytobands, URL regions) throws FileNotFoundException, URISyntaxException {
		super(cytobands, CytobandHandlerThread.class);
		this.regionDataSource = new LineDataSource(regions, CytobandHandlerThread.class);
	}
	
	public LineDataSource getRegionDataSrouce() {
		return regionDataSource;
	}

}
