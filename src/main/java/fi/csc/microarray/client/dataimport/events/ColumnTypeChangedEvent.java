package fi.csc.microarray.client.dataimport.events;

import fi.csc.microarray.client.dataimport.ColumnType;

public class ColumnTypeChangedEvent extends ImportChangeEvent {

	private int columnIndex;
	
	public ColumnTypeChangedEvent(Object source, ColumnType newValue, int columnIndex) {
		super(source, newValue);
		this.columnIndex = columnIndex;
	}
	
	public int getColumnIndex(){
		return this.columnIndex;
	}

	@Override
	public ColumnType getNewValue() {
		return (ColumnType)this.newValue;
	}
}
