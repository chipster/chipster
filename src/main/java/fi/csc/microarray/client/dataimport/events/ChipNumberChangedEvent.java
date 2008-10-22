package fi.csc.microarray.client.dataimport.events;


public class ChipNumberChangedEvent extends ImportChangeEvent {

	private int columnIndex;
	
	public ChipNumberChangedEvent(Object source, Integer newValue, int columnIndex) {
		super(source, newValue);
		this.columnIndex = columnIndex;
	}

	@Override
	public Integer getNewValue() {
		return (Integer)this.newValue;
	}

	public int getColumnIndex(){
		return this.columnIndex;
	}
	
}
