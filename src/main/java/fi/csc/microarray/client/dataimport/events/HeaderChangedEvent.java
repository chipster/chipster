package fi.csc.microarray.client.dataimport.events;


public class HeaderChangedEvent extends TableRowChangeEvent {

	public HeaderChangedEvent(Object source, int newValue) {
		super(source, newValue);
	}
}
