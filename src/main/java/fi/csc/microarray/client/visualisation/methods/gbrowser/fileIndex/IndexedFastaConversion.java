package fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex;

import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;

public class IndexedFastaConversion extends DataThread {

	private IndexedFastaDataSource dataSource;

	public IndexedFastaConversion(URL data, URL index, final GBrowser browser) {

		super(browser);

		try {			
			this.dataSource = new IndexedFastaDataSource(data, index);

		} catch (URISyntaxException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Override
	public void clean() {
	}


	@Override
	protected void processDataRequest(DataRequest request) {
		
		if (request.start.bp < 1) {
			request.start.bp = 1l;
		}

		List<RegionContent> responseList = new LinkedList<RegionContent>();
		
		String sequence = dataSource.query(request.start.chr, request.start.bp, request.end.bp);

		LinkedHashMap<DataType, Object> values = new LinkedHashMap<DataType, Object>();
		values.put(DataType.SEQUENCE, sequence);

		RegionContent regCont = new RegionContent(request, values);

		/*
		 * NOTE! RegionContents created from the same read area has to be equal in methods equals, hash and compareTo. Primary types
		 * should be ok, but objects (including tables) has to be handled in those methods separately. Otherwise tracks keep adding
		 * the same reads to their read sets again and again.
		 */
		responseList.add(regCont);

		// Send result
		createDataResult(new DataResult(request.getStatus(), responseList));				
	}

	public String toString() {
		return this.getClass().getName() + " - " + dataSource;
	}
}
