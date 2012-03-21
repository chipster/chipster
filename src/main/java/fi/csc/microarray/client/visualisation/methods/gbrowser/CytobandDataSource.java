package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.FileNotFoundException;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.CytobandHandlerThread;

public class CytobandDataSource extends LineDataSource {
	
	private LineDataSource regionDataSource;
	private LineDataSource coordDataSource;

	public CytobandDataSource(URL cytobands, URL regions, URL coordSystem) throws FileNotFoundException, URISyntaxException {
		super(cytobands, CytobandHandlerThread.class);
		this.regionDataSource = new LineDataSource(regions, CytobandHandlerThread.class);
		this.coordDataSource = new LineDataSource(coordSystem, CytobandHandlerThread.class);
	}
	
	public LineDataSource getRegionDataSource() {
		return regionDataSource;
	}
	
	public LineDataSource getCoordDataSource() {
		return coordDataSource;
	}
}
