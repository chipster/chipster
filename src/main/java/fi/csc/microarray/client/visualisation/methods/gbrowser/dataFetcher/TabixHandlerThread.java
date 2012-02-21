package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.IOException;
import java.util.Queue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;

/**
 * Experimental Tabix support.
 * 
 * @author Aleksi Kallio
 *
 */
public class TabixHandlerThread extends AreaRequestHandler {

	TabixDataSource tabixData;

	public TabixHandlerThread(DataSource file, Queue<AreaRequest> areaRequestQueue,
			AreaResultListener areaResultListener) {

		super(areaRequestQueue, areaResultListener);
		tabixData = (TabixDataSource) file;
	}

	/**
	 * Handles normal and concised area requests by using TabixFile.
	 */
	@Override
	protected void processAreaRequest(AreaRequest areaRequest) {       
		try {
			createAreaResult(new AreaResult(areaRequest.status, tabixData.getTabix().getReads(areaRequest)));
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
