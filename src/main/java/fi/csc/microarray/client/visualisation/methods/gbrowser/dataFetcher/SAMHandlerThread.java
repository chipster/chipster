package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.Queue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.SAMDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class SAMHandlerThread extends AreaRequestHandler {
    
    SAMDataSource samData;

    public SAMHandlerThread(DataSource file, Queue<AreaRequest> areaRequestQueue,
            AreaResultListener areaResultListener) {
        
        super(areaRequestQueue, areaResultListener);
        samData = (SAMDataSource) file;
    }

    /**
     * Handles normal and concised area requests by using SAMFile.
     */
    @Override
    protected void processAreaRequest(AreaRequest areaRequest) {       
        if (areaRequest.status.concise) {
            // Create concise results
            for (RegionContent content : samData.getSAM().getConciseReads(areaRequest)) {
                createAreaResult(new AreaResult<RegionContent>(areaRequest.status, content));
            }
        } else {
            // Create a result for each read
            for (RegionContent content : samData.getSAM().getReads(areaRequest)) {
                createAreaResult(new AreaResult<RegionContent>(areaRequest.status, content));
            }            
        }
    }

}
