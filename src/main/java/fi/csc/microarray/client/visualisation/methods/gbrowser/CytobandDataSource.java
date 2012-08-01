package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.FileNotFoundException;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.CytobandHandlerThread;

public class CytobandDataSource extends LineDataSource {

	public CytobandDataSource(URL cytobands) throws FileNotFoundException, URISyntaxException {
		super(cytobands, CytobandHandlerThread.class);
	}
}
