package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.RegionContent;


public interface AreaResultListener {
	public void processAreaResult(AreaResult<RegionContent> areaResult);
}
