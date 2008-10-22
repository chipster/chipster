package fi.csc.microarray.client.dataimport.events;

import fi.csc.microarray.client.dataimport.Delimiter;


public class DelimiterChangedEvent extends ImportChangeEvent {

	public DelimiterChangedEvent(Object source, Delimiter newValue) {
		super(source, newValue);
	}

	public Delimiter getNewValue(){
		return (Delimiter)this.newValue;
	}
}
