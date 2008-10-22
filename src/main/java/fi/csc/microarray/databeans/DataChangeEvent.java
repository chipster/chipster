package fi.csc.microarray.databeans;

public abstract class DataChangeEvent {
    	
	private DataItem dataItem;
	
	public DataChangeEvent(DataItem dataItem) {
		this.dataItem = dataItem;
	}

	public DataItem getDataItem() {
		return dataItem;
	}
}
