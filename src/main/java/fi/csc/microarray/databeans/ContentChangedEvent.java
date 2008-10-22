package fi.csc.microarray.databeans;

public class ContentChangedEvent extends DataChangeEvent {

	public ContentChangedEvent(DataItem chosenOne) {
		super(chosenOne);
	}
}
