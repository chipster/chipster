package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.Queue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.SAMDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class SAMHandlerThread extends AreaRequestHandler {
    
    SAMDataSource samData;

    public SAMHandlerThread(DataSource file, Queue<AreaRequest> areaRequestQueue,
            AreaResultListener areaResultListener) {
        
        super(areaRequestQueue, areaResultListener);
        samData = (SAMDataSource) file;
    }

    @Override
    protected void processAreaRequest(AreaRequest areaRequest) {     
        FsfStatus status = areaRequest.status;
        
        // Create a result for each read
        for (RegionContent content : samData.getSAM().getReads(areaRequest)) {
            createAreaResult(new AreaResult<RegionContent>(status, content));
        }
    }

}
