package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class GeneSearchConversion extends SingleThreadAreaRequestHandler {	
	
	private HashMap<String, Chromosome> geneNameMap;
	private LineDataSource dataSource;
	
	public GeneSearchConversion(URL url, final GBrowser browser) {

		super(null, null);

		try {

			dataSource = new LineDataSource(url, null);

		} catch (IOException e) {
			e.printStackTrace();
		} catch (URISyntaxException e) {
			e.printStackTrace();
		}
	}

	@Override
	protected void processAreaRequest(AreaRequest request) {

		super.processAreaRequest(request);	

		if (request.getStatus().poison) {
			return;
		}

		if (geneNameMap == null) {
			geneNameMap = new HashMap<String, Chromosome>();

			try {
				readFile();
				
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		GeneRequest geneRequest = (GeneRequest)request;
		
		String searchString = geneRequest.getSearchString();

		Chromosome chr = geneNameMap.get(searchString.toLowerCase());
		List<RegionContent> resultList = new LinkedList<RegionContent>();

		if (chr != null) {

			resultList.add(new RegionContent(new Region(null, null, chr), null));
		}

		super.createAreaResult(new GeneResult(geneRequest.getStatus(), resultList, searchString));
	}
	
	private void readFile() throws IOException {

		String line;

		while ((line = dataSource.readLine()) != null) {

			String[] cols = line.split("\t");

			if (cols.length == 2) {
				String chr = cols[0];
				String geneName = cols[1];

				geneNameMap.put(geneName.toLowerCase(), new Chromosome(chr));
			}
		}
	}
}
