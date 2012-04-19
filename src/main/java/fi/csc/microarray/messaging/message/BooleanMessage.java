package fi.csc.microarray.messaging.message;

import javax.jms.JMSException;
import javax.jms.MapMessage;


public class BooleanMessage extends ChipsterMessage {
	
	private final static String KEY_BOOLEAN = "boolean";
	
	private boolean value;
	
	public BooleanMessage() {
		super();
	}
	
	public BooleanMessage(boolean value) {
		super();
		this.value = value;
	}
	
	public boolean getValue() {
		return value;
	}
	
	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);
		this.value = Boolean.parseBoolean(from.getString(KEY_BOOLEAN));
	}

	public void marshal(MapMessage mapMessage) throws JMSException {
		super.marshal(mapMessage);
		mapMessage.setString(KEY_BOOLEAN, String.valueOf(this.value));
	}
}
