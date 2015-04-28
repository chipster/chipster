/*
 * Created on Feb 11, 2005
 *
 *
 */
package fi.csc.microarray.messaging.message;

import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.jms.JMSException;
import javax.jms.MapMessage;

import org.apache.log4j.Logger;


/**
 * For sending messages with payloads. Payload means any data whose structure
 * is not defined in the JMSMessage. For example data bean contents and tool
 * descriptions are transferred as payloads.
 * 
 * The actual data is sent through the FileBroker and payload message only
 * carries the id for the data in the FileBroker.
 * 
 * @author Taavi Hupponen, Aleksi Kallio
 *
 */
public class PayloadMessage extends ParameterMessage {

	private static final Logger logger = Logger.getLogger(PayloadMessage.class);

	private static final String KEY_ID_PREFIX = "payload_id_";
	private static final String KEY_NAME_PREFIX = "payload_name_";
	
	private Map<String, String> ids = new HashMap<String, String>();
	private Map<String, String> names = new HashMap<String, String>();
	
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

		// load payload ids
		try {
			for (Enumeration<String> keys = from.getMapNames(); keys.hasMoreElements(); ) {
				String name = keys.nextElement();
				logger.debug("examining " + name);
				if (name.startsWith(KEY_ID_PREFIX)) {
					String payloadName = name.substring(KEY_ID_PREFIX.length());
					String value  = from.getString(name);
					ids.put(payloadName, value);
					logger.debug("Unmarshalled " + name + " -> " + payloadName + ", " + value);
				}
				if (name.startsWith(KEY_NAME_PREFIX)) {
					String payloadName = name.substring(KEY_NAME_PREFIX.length());
					String value  = from.getString(name);
					names.put(payloadName, value);
					logger.debug("Unmarshalled " + name + " -> " + payloadName + ", " + value);
				}
			}
		} catch (Exception e) {
			handleException(e);
		}
	}

	
	@Override
	public void marshal(MapMessage mapMessage) throws JMSException {
		super.marshal(mapMessage);
		
		// add payload ids
		String key;
		String idString;
		try {
			for (String name : ids.keySet()) {
				key = KEY_ID_PREFIX + name;
				idString = ids.get(name);
				mapMessage.setString(key, idString);
				logger.debug("Marshalled " + name + " -> " + key + " " + idString);
			}
			for (String name : names.keySet()) {
				key = KEY_NAME_PREFIX + name;
				idString = names.get(name);
				mapMessage.setString(key, idString);
				logger.debug("Marshalled " + name + " -> " + key + " " + idString);
			}
		} catch (Exception e) {
			handleException(e);
		}
		
	}

	
	/**
	 * Add the id of an already uploaded payload to the payloads.
	 * 
	 * 
	 * @param key name of the payload (the input name of the operation, not the name of the databean)
	 * @param id the id of the payload on the filebroker
	 * @param name the dataset name
	 */
	public void addPayload(String key, String id, String name) {
		ids.put(key, id);
		names.put(key, name);
	}


	/**
	 * Return the id for the payload content.
	 * 
	 * @param name
	 * @return
	 * @throws JMSException
	 */
	public String getId(String key) throws JMSException {
		if (ids.containsKey(key)) {
			return ids.get(key);
		} else {
			throw new IllegalArgumentException("No payload with name: " + key);
		}
	}

	/**
	 * Return the name for the payload content.
	 * 
	 * @param name
	 * @return
	 * @throws JMSException
	 */
	public String getName(String key) throws JMSException {
		if (names.containsKey(key)) {
			return names.get(key);
		} else {
			throw new IllegalArgumentException("No payload with name: " + key);
		}
	}
	
	/**
	 * Get all payload keys.
	 * 
	 * @return
	 */
	public Set<String> getKeys() {
		HashSet<String> keys = new HashSet<>();
		keys.addAll(ids.keySet());
		keys.addAll(names.keySet());
		return keys;
	}
	
	@Override
	public String toString() {
		return super.toString() + ", payload count: " + getKeys().size();		
	}

}
