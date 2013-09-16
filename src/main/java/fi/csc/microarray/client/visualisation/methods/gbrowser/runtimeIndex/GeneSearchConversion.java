package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

/**
 * This class converts gene name search requests to GeneResult objects containing 
 * the chromosome of the requested gene. This information is read from the purpose-build
 * file. 
 * 
 * @author klemela
 */
public class GeneSearchConversion extends DataThread {	
	
	private HashMap<String, Chromosome> geneNameMap;
	private LineDataSource dataSource;
	
	public GeneSearchConversion(DataUrl data, final GBrowser browser) {

		super(browser, null);

		try {

			dataSource = new LineDataSource(data);
			super.setDataSource(dataSource);

		} catch (IOException e) {
			e.printStackTrace();
		} catch (URISyntaxException e) {
			e.printStackTrace();
		}
	}

	@Override
	protected void processDataRequest(DataRequest request) {

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
		List<Feature> resultList = new LinkedList<Feature>();

		if (chr != null) {

			resultList.add(new Feature(new Region(null, null, chr), null));
		}

		super.createDataResult(new GeneResult(geneRequest.getStatus(), resultList, searchString));
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
