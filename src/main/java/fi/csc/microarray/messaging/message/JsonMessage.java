package fi.csc.microarray.messaging.message;

import javax.jms.JMSException;
import javax.jms.MapMessage;

public class JsonMessage extends ChipsterMessage {
	
	private final static String KEY_JSON = "json";
	
	private String value;
	
	public JsonMessage() {
		super();
	}
	
	public JsonMessage(String value) {
		super();
		this.value = value;
	}
	
	public String getJson() {
		return value;
	}
	
	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);
		this.value = from.getString(KEY_JSON);
	}

	public void marshal(MapMessage mapMessage) throws JMSException {
		super.marshal(mapMessage);
		mapMessage.setString(KEY_JSON, this.value);
	}
}