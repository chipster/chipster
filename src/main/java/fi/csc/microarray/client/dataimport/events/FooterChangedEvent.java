package fi.csc.microarray.client.dataimport.events;


public class FooterChangedEvent extends TableRowChangeEvent {

	public FooterChangedEvent(Object source, int newValue) {
		super(source, newValue);
	}
}
