package fi.csc.microarray.client.dataimport.events;


public class ChipCountChangeEvent extends ImportChangeEvent {
	
	public ChipCountChangeEvent(Object source, Integer newValue) {
		super(source, newValue);
	}

	@Override
	public Integer getNewValue() {
		return (Integer)this.newValue;
	}
}
