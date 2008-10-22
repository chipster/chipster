package fi.csc.microarray.client.dataimport.events;


public class TableRowChangeEvent extends ImportChangeEvent {

	public TableRowChangeEvent(Object source, Object newValue) {
		super(source, newValue);
	}

	public Integer getNewValue(){
		return (Integer)this.newValue;
	}

}
