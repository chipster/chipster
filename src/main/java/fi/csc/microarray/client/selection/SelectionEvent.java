package fi.csc.microarray.client.selection;

import java.beans.PropertyChangeEvent;

import fi.csc.microarray.databeans.DataBean;



/**
 * <p>An event that is fired whenever selection changes. The event
 * contains only information that is needed in deciding if
 * the component is interested in this event. To find out
 * what the current selection is, the component should 
 * access IntegratedSelectionManager.</p>
 * 
 * <p>SelectionEvents are sent by IntegratedSelectionManager.
 * It tracks sources so that the source of the
 * event is not the manager, but the component that actually
 * triggered the selection change.</p> 
 * 
 * @see IntegratedSelectionManager
 * 
 * @author Petri Klemelä, Aleksi Kallio
 *
 */
public class SelectionEvent extends PropertyChangeEvent {
	private DataBean data;
	
	public SelectionEvent(DataBean data, Object source) {
		super(source, null, null, null);
		this.data = data;
	}
	
	public DataBean getData(){
		return data;
	}
}
