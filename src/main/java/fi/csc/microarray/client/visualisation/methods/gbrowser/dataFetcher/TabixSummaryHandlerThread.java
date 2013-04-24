package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.IOException;
import java.util.Queue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;

/**
 * Experimental Tabix support.
 * 
 * @author Aleksi Kallio
 *
 */
public class TabixSummaryHandlerThread extends AreaRequestHandler {

	TabixSummaryDataSource tabixData;

	public TabixSummaryHandlerThread(DataSource file, Queue<AreaRequest> areaRequestQueue,
			AreaResultListener areaResultListener) {

		super(areaRequestQueue, areaResultListener);
		tabixData = (TabixSummaryDataSource) file;
	}

	public TabixSummaryHandlerThread(DataSource file) {
		this(file, null, null);
	}

	/**
	 * Handles normal and concised area requests by using TabixFile.
	 */
	@Override
	protected void processAreaRequest(AreaRequest areaRequest) {       
		
		//No other threads to poison, this one killed in parent class.
		super.processAreaRequest(areaRequest);
		
		if (areaRequest.getStatus().poison) {
			return;
		}
		
		try {
			createAreaResult(new AreaResult(areaRequest.getStatus(), tabixData.getTabix().getReads(areaRequest)));
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
