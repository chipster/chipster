package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.File;
import java.io.FileNotFoundException;

public class CytobandDataSource extends LineDataSource {
	
	private LineDataSource regionDataSource;

	public CytobandDataSource(File cytobandFile, File regionFile) throws FileNotFoundException {
		super(cytobandFile);
		this.regionDataSource = new LineDataSource(regionFile);
	}
	
	public LineDataSource getRegionDataSrouce() {
		return regionDataSource;
	}

}
