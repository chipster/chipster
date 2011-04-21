package fi.csc.microarray.client.session;

import java.util.LinkedList;
import java.util.List;

import javax.xml.bind.ValidationEvent;
import javax.xml.bind.ValidationEventHandler;

public class NonStoppingValidationEventHandler implements ValidationEventHandler {
	
	private List<ValidationEvent> validationEvents = new LinkedList<ValidationEvent>();
	
	/**
	 * Continue, no matter what.
	 */
	@Override
	public boolean handleEvent(ValidationEvent event) {
		this.validationEvents.add(event);
		return true;
	}
	
	public boolean hasEvents() {
		return validationEvents.size() > 0;
	}

	public String getValidationEventsAsString() {
		String s = "";
		for (ValidationEvent event : validationEvents) {
			s += event.getMessage() + "\n";
		}
		return s;
	}
}
