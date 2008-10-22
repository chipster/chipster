package fi.csc.microarray.databeans;

public class DataItemRemovedEvent extends DataChangeEvent {

	public DataItemRemovedEvent(DataItem dataItem) {
		super(dataItem);
	}
}
