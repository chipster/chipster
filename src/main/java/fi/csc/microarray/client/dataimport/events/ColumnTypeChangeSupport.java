package fi.csc.microarray.client.dataimport.events;

public interface ColumnTypeChangeSupport {

	public void addColumnTypeChangeListener(ColumnTypeChangeListener l);
	public void removeColumnTypeChangeListener(ColumnTypeChangeListener l);
	public void fireColumnTypeChangeEvent(ColumnTypeChangedEvent event);
	public void fireChipNumberChangeEvent(ChipNumberChangedEvent event);
	public void fireChipCountChangeEvent(ChipCountChangeEvent event);
	
}
