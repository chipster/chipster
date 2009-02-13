package fi.csc.microarray.messaging.message;

import java.net.MalformedURLException;
import java.net.URL;

import javax.jms.JMSException;
import javax.jms.MapMessage;


public class UrlMessage extends NamiMessage {
	
	private final static String KEY_URL = "url";
	
	private URL url;
	
	public UrlMessage(URL url) {
		super();
		this.url = url;
	}
	
	public URL getUrl() {
		return url;
	}
	
	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);
		try {
			this.url = new URL(from.getString(KEY_URL));
		} catch (MalformedURLException e) {
			handleException(e);
		}
	}

	public void marshal(MapMessage mapMessage) throws JMSException {
		super.marshal(mapMessage);
	}
}
