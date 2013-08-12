package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;

import org.broad.tribble.readers.TabixReader;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class BedTabixToRegionConversion extends DataThread {

	private TabixDataSource dataSource;
	private BedLineParser parser = new BedLineParser(false);

	public BedTabixToRegionConversion(URL tabixFile, URL tabixIndexFile, final GBrowser browser) {

		super(browser);

		try {
			this.dataSource = new TabixDataSource(tabixFile, tabixIndexFile);

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
			List<RegionContent> resultList = new LinkedList<RegionContent>();

			if (iterator != null) { //null if there isn't such chromosome in annotations

				while ((line = iterator.next()) != null) {

					parser.setLine(line);
					Region region = parser.getRegion();				
					resultList.add(new RegionContent(region));
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
