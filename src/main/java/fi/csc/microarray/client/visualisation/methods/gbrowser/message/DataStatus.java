package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;

/**
 * <p>Data retrieval status</p>
 * 
 * <p>A generic status field that travels with all requests.</p> 
 */
public class DataStatus implements Cloneable {
	
	public DataStatus() {		
	}
	
	public DataStatus(DataStatus dataStatus) {
		this.dataRequestCount = dataStatus.dataRequestCount;
		this.dataThread = dataStatus.dataThread;
	}

	private long dataRequestCount = -1;
	private DataThread dataThread;
	
	public DataThread getDataThread() {
		return dataThread;
	}
	public void setDataThread(DataThread dataThread) {
		this.dataThread = dataThread;
	}
	public long getDataRequestCount() {
		return dataRequestCount;
	}
	public void setDataRequestCount(long dataRequestCount) {
		this.dataRequestCount = dataRequestCount;
	}
}
