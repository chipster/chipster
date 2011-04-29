package fi.csc.microarray.messaging.message;

import javax.jms.JMSException;
import javax.jms.MapMessage;

public class SourceMessage extends ChipsterMessage {
	
	private final static String KEY_SOURCE = "source";
	
	private String source;
	
	public SourceMessage() {
		super();
	}
	
	public SourceMessage(String source) {
		super();
		this.source = source;
	}
	
	public String getSource() {
		return source;
	}
	
	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);
			this.source = from.getString(KEY_SOURCE);
	}

	public void marshal(MapMessage mapMessage) throws JMSException {
		super.marshal(mapMessage);
		mapMessage.setString(KEY_SOURCE, source);
	}
}
