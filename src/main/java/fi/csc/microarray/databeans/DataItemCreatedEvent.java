package fi.csc.microarray.databeans;

public class DataItemCreatedEvent extends DataChangeEvent {

	public DataItemCreatedEvent(DataItem dataItem) {
		super(dataItem);
	}
}
