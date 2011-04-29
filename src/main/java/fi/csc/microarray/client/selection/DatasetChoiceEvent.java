package fi.csc.microarray.client.selection;

import java.beans.PropertyChangeEvent;



/**
 * An event that is fires whenever a DataBean has been chosen in a
 * DataSetView component. Information about the selection will then be
 * forwarded to other DataSetViews so that they can respond accordingly.
 * 
 * @author Janne Käki, Aleksi Kallio
 *
 */
public class DatasetChoiceEvent extends PropertyChangeEvent {

	public DatasetChoiceEvent(Object source) {
		super(source, null, null, null);
	}
}
