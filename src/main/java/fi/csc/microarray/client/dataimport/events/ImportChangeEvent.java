package fi.csc.microarray.client.dataimport.events;

import java.util.EventObject;

/**
 * Abstract root class from which all import event should be derived
 * 
 * @author mkoski
 *
 */
public abstract class ImportChangeEvent extends EventObject {

	protected Object newValue;
	
	public ImportChangeEvent(Object source, Object newValue) {
		super(source);
		this.newValue = newValue;
	}

	public abstract Object getNewValue();
	
	public String toString(){
		return this.getClass().getName() + " newValue: " + this.newValue.toString();
	}
}
