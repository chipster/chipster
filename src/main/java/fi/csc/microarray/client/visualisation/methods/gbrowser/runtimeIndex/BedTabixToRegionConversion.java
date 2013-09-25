package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.LinkedList;
import java.util.List;

import org.broad.tribble.readers.TabixReader;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class BedTabixToRegionConversion extends DataThread {

	private TabixDataSource dataSource;
	private BedLineParser parser = new BedLineParser(false);

	public BedTabixToRegionConversion(DataUrl repeatUrl, DataUrl repeatIndexUrl, final GBrowser browser) {

		super(browser, null);

		try {
			this.dataSource = new TabixDataSource(repeatUrl, repeatIndexUrl);
			super.setDataSource(dataSource);

		} catch (URISyntaxException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Override
	public void clean() {

		dataSource.clean();
	}


	@Override
	protected void processDataRequest(DataRequest request) {

		// Read the given region
		TabixReader.Iterator iterator = dataSource.getTabixIterator(request);

		try {
			String line;
			List<Feature> resultList = new LinkedList<Feature>();

			if (iterator != null) { //null if there isn't such chromosome in annotations

				while ((line = iterator.next()) != null) {

					parser.setLine(line);
					Region region = parser.getRegion();				
					resultList.add(new Feature(region));
				}
			}

			// Send result			
			super.createDataResult(new DataResult(request.getStatus(), resultList));		
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public String toString() {
		return this.getClass().getName() + " - " + dataSource;
	}
}
