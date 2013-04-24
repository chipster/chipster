package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;

import org.broad.tribble.readers.TabixReader;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.SamBamUtils;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.TabixUtil;

public class BedTabixToRegionConversion extends SingleThreadAreaRequestHandler {

	private TabixDataSource dataSource;
	private BedLineParser parser = new BedLineParser(false);

	public BedTabixToRegionConversion(URL tabixFile, URL tabixIndexFile, final GBrowser browser) {

		super(null, null);

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

		SamBamUtils.closeIfPossible(dataSource.getReader());
	}


	@Override
	protected void processAreaRequest(AreaRequest request) {

		super.processAreaRequest(request);

		if (request.getStatus().poison) {
			return;
		}


		// Read the given region
		TabixReader.Iterator iterator = TabixUtil.getTabixIterator(dataSource, request);

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
			super.createAreaResult(new AreaResult(request.getStatus(), resultList));		
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public String toString() {
		return this.getClass().getName() + " - " + dataSource;
	}
}
