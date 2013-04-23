package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

public class CytobandConversion extends SingleThreadAreaRequestHandler {

	private Index index;

	private CytobandLineParser parser;

	public CytobandConversion(DataSource file, final GBrowser browser) {

		super(null, null);

		this.parser = new CytobandLineParser();

		try {

			this.index = new InMemoryIndex(file, parser);

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Override
	protected void processAreaRequest(AreaRequest request) {

		super.processAreaRequest(request);	

		if (request.getStatus().poison) {
			return;
		}

		if (index == null) {
			return;
		}


		TreeMap<IndexKey, String> lineMap;
		try {
			lineMap = index.getFileLines(request);

			List<RegionContent> resultList = new ArrayList<RegionContent>();

			for (Entry<IndexKey, String> entry : lineMap.entrySet()) {

				LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();

				parser.setLine(entry.getValue());				
				Region region = parser.getRegion();
				Cytoband cytoband = new Cytoband(region, parser.getBand(), parser.getStain());

				values.put(ColumnType.ID, entry.getKey());
				values.put(ColumnType.VALUE, cytoband);

				resultList.add(new RegionContent(region, values));
			}

			super.createAreaResult(new AreaResult(request.getStatus(), resultList));

		} catch (IOException e) {
			e.printStackTrace();
		} catch (GBrowserException e) {
			e.printStackTrace();
		}
	}
}
