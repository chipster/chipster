package fi.csc.microarray.client.dataimport.events;

/**
 * Event for column title value change. Note that there is also own method 
 * for title row number change.
 * 
 * @author mkoski
 *
 */
public class ColumnTitlesChangedEvent extends ImportChangeEvent {

	public ColumnTitlesChangedEvent(Object source, String[] newValue) {
		super(source, newValue);
	}

	@Override
	public String[] getNewValue() {
		return (String[])newValue;
	}

}
