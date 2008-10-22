package fi.csc.microarray.client.dataimport.events;


public class DecimalSeparatorChangedEvent extends ImportChangeEvent {
	
	public DecimalSeparatorChangedEvent(Object source, Character newValue) {
		super(source, newValue);
	}
	
	public Character getNewValue(){
		return (Character)this.newValue;
	}

}
