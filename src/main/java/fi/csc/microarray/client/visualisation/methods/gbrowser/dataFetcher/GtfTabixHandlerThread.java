package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.Queue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ParsedFileResult;


public class GtfTabixHandlerThread extends TabixHandlerThread {
    

	private GtfTabixFileFetcherThread fileFetcher;

	
	public GtfTabixHandlerThread(DataSource file,
			Queue<AreaRequest> areaRequestQueue,
			AreaResultListener areaResultListener) {
		super(file, areaRequestQueue, areaResultListener);
	}

	public GtfTabixHandlerThread(TabixDataSource gtfDataSource) {
		this(gtfDataSource, null, null);
	}

	@Override
	public void runThread() {

		// Start file processing layer thread
		fileFetcher = new GtfTabixFileFetcherThread(fileRequestQueue, fileResultQueue, this, dataSource);
		fileFetcher.start();
		
		// Start this thread
		super.runThread();
	}

	@Override
	protected void processFileResult(ParsedFileResult fileResult) {

		if (fileResult.getFileRequest().areaRequest instanceof GeneRequest) {
			GeneRequest geneRequest = (GeneRequest)fileResult.getFileRequest().areaRequest;
			createAreaResult(new GeneResult(fileResult.getStatus(), fileResult.getContents(), geneRequest.getSearchString()));
		} else {
			createAreaResult(new AreaResult(fileResult.getStatus(), fileResult.getContents()));
		}
	}
}
