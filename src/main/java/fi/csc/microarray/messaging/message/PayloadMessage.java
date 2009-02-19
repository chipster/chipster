/*
 * Created on Feb 11, 2005
 *
 *
 */
package fi.csc.microarray.messaging.message;

import java.net.URL;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.jms.JMSException;
import javax.jms.MapMessage;

import org.apache.log4j.Logger;


/**
 * For sending jobs to back-end components.
 * 
 * @author Taavi Hupponen, Aleksi Kallio
 *
 */
public class PayloadMessage extends ParameterMessage {

	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(PayloadMessage.class);

	private static final String KEY_PAYLOAD_PREFIX = "payload_";
	
	private Map<String, URL> payloads = new HashMap<String, URL>();
	
	/**
	 * For reflection compatibility (newInstance). DO NOT REMOVE!
	 */
	public PayloadMessage() {
		super();
	}
	
	public PayloadMessage(List<String> parameters) {
		super(parameters);
	}

	@SuppressWarnings("unchecked")
	@Override
	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);

		// load payload urls
		try {
			for (Enumeration<String> names = from.getMapNames(); names.hasMoreElements(); ) {
				String name = names.nextElement();
				logger.debug("examining " + name);
				if (name.startsWith(KEY_PAYLOAD_PREFIX)) {
					String payloadName = name.substring(KEY_PAYLOAD_PREFIX.length());
					URL url = new URL(from.getString(name));
					payloads.put(payloadName, url);
					logger.debug("Unmarshalled " + name + " -> " + payloadName + ", " + url.toExternalForm());
				}
			}
		} catch (Exception e) {
			handleException(e);
		}
	}

	
	/**
	 * Construct a MapMessage that can be used to create a new JobMessage.
	 */
	@Override
	public void marshal(MapMessage mapMessage) throws JMSException {
		super.marshal(mapMessage);
		
		// add payload urls
		String key;
		String urlString;
		try {
			for (String name : payloadNames()) {
				key = KEY_PAYLOAD_PREFIX + name;
				urlString = payloads.get(name).toExternalForm();
				mapMessage.setString(key, urlString);
				logger.debug("Marshalled " + name + " -> " + key + " " + urlString);
			}
		} catch (Exception e) {
			handleException(e);
		}
		
	}

	
	/**
	 * Add the URL of an already uploaded payload to the payloads.
	 * 
	 * 
	 * @param payloadName name of the payload (the input name of the operation, not the name of the databean)
	 * @param payloadURL an URL pointing to the server side location of the payload
	 */
	public void addPayload(String payloadName, URL payloadURL) {
		payloads.put(payloadName, payloadURL);
	}


	/**
	 * Return the URL for the payload content.
	 * 
	 * @param name
	 * @return
	 * @throws JMSException
	 */
	public URL getPayload(String name) throws JMSException {
		if (payloads.containsKey(name)) {
			return payloads.get(name);
		} else {
			throw new IllegalArgumentException("No payload with name: + name");
		}
	}
	
	
	public Set<String> payloadNames() {
		return payloads.keySet();
	}


	public Map<String, URL> getPayloads() {
		return payloads;
	}

	
	public void setPayloads(Map<String, URL> payloads) {
		this.payloads = payloads;
	}

	
	@Override
	public String toString() {
		return super.toString() + ", payload count: " + payloadNames().size();		
	}
	


}
