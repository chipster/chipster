package fi.csc.microarray.client.selection;

import java.beans.PropertyChangeEvent;

import fi.csc.microarray.databeans.Dataset;



/**
 * An event that is fired whenever data rows are selected from data set 
 * to keep the visualisations agreed on the selected rows. 
 * 
 * @author Petri KlemelÃ¤
 *
 */
public class RowChoiceEvent extends PropertyChangeEvent {
	Dataset data;
	
	public RowChoiceEvent(Dataset data, Object source) {
		super(source, null, null, null);
		this.data = data;
	}
	
	public Dataset getData(){
		return data;
	}
}
