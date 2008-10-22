package fi.csc.microarray.client.selection;

import java.beans.PropertyChangeEvent;

import fi.csc.microarray.databeans.DataBean;



/**
 * An event that is fired whenever data rows are selected from data set 
 * to keep the visualisations agreed on the selected rows. 
 * 
 * @author Petri Klemelä
 *
 */
public class RowChoiceEvent extends PropertyChangeEvent {
	DataBean data;
	
	public RowChoiceEvent(DataBean data, Object source) {
		super(source, null, null, null);
		this.data = data;
	}
	
	public DataBean getData(){
		return data;
	}
}
