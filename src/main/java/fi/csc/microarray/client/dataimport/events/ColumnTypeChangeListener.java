package fi.csc.microarray.client.dataimport.events;

import java.util.EventListener;

/**
 * Listener for column type changes and chip number changes
 * @author mkoski
 *
 */
public interface ColumnTypeChangeListener extends EventListener {

	public void columnTypeChanged(ColumnTypeChangedEvent event);
	public void chipNumberChanged(ChipNumberChangedEvent event);
	public void chipCountChanged(ChipCountChangeEvent event);
}
