package fi.csc.microarray.client.dataimport.events;

import java.io.File;

public class InputFileChangedEvent extends ImportChangeEvent {

	public InputFileChangedEvent(Object source, File newValue) {
		super(source, newValue);
	}

	@Override
	public File getNewValue() {
		return (File)getNewValue();
	}

}
