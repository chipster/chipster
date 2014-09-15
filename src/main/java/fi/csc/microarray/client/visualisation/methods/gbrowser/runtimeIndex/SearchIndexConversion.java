package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.SearchRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

/**
 * Search indexes are TSV files with columns chr, start, end and key. Search 
 * finds the row where the key equals search term (ignoring case) and returns 
 * its region. An instance is initialized with a list of URLs. When there are 
 * multiple URLs, those are searched in the list's order until a matching key 
 * is found.
 * 
 * At the moment, there are indexes only for gene names and ids, These files 
 * are small, only a few megabytes and are simply read through and saved in 
 * memory.
 *  
 * To prepare for the future needs, search indexes are sorted according to 
 * key (ignoring case). It is possible to implement binary search for keys 
 * similarly to how genome browser usually locates positions in TSV, GTF and 
 * BED files, if larger index files are needed at some point.
 * 
 * @author klemela
 */
public class SearchIndexConversion extends DataThread {
	
	public static class SearchIndex {
		
		public SearchIndex(DataUrl dataUrl) throws URISyntaxException, IOException {
			dataSource = new LineDataSource(dataUrl);
		}
		
		private LineDataSource dataSource;
		private HashMap<String, Region> indexMap;
		
		public Region search(String searchString) {
			if (indexMap == null) {				
				indexMap = new HashMap<String, Region>();

				try {
					readFile();				
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			
			return indexMap.get(searchString);
		}
		
		private void readFile() throws IOException {

			String line;

			while ((line = dataSource.readLine()) != null) {

				String[] cols = line.split("\t");

				String chr = cols[0];
				long start = Long.parseLong(cols[1]);
				long end = Long.parseLong(cols[2]);
				String string = cols[3];

				indexMap.put(string.toLowerCase(), new Region(start, end, new Chromosome(chr)));			
			}
		}
	}
	
	ArrayList<SearchIndex> indexes = new ArrayList<>();
	
	public SearchIndexConversion(List<DataUrl> indexUrls, final GBrowser browser) {

		super(browser, null);

		try {

			for (DataUrl dataUrl : indexUrls) {
				indexes.add(new SearchIndex(dataUrl));
			}			

		} catch (IOException e) {
			e.printStackTrace();
		} catch (URISyntaxException e) {
			e.printStackTrace();
		}
	}

	@Override
	protected void processDataRequest(DataRequest req) throws InterruptedException {

		SearchRequest request = (SearchRequest)req;
		
		String searchString = request.getSearchString().toLowerCase();
		
		for (SearchIndex index : indexes) {
			Region region = index.search(searchString);
			if (region != null) {				
				reply(request, region);											
				return;
			}
		}		
		// not found
		reply(request, null);
	}

	private void reply(SearchRequest request, Region region) throws InterruptedException {
		List<Feature> resultList = new LinkedList<Feature>();							
		resultList.add(new Feature(region, null));
		super.createDataResult(new GeneResult(
				request.getStatus(), resultList, request.getSearchString()));
	}
}
