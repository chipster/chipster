package fi.csc.microarray.client.dataimport.events;

/**
 * Event for title header row change. Notice that this event means 
 * that the title header <strong>row number</strong> is changed. There 
 * is also an event for column title values change.
 * 
 * @author mkoski
 *
 */
public class TitleRowChangedEvent extends TableRowChangeEvent {
	
	public TitleRowChangedEvent(Object source, int newValue) {
		super(source, newValue);
	}
}
