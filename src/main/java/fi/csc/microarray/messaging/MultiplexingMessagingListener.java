package fi.csc.microarray.messaging;

import org.apache.log4j.Logger;

import java.util.HashMap;
import java.util.Map;

import fi.csc.microarray.messaging.message.ChipsterMessage;

public class MultiplexingMessagingListener implements MessagingListener {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger
			.getLogger(MultiplexingMessagingListener.class);

	private Map<String, MessagingListener> channels = new HashMap<String, MessagingListener>();
	
	public void addChannel(String channelName, MessagingListener channelListener) {
		channels.put(channelName, channelListener);
	}
	
	public void removeChannel(String channelName) {
		channels.remove(channelName);
	}
	
	public void onChipsterMessage(ChipsterMessage msg) {
		String channelName = msg.getMultiplexChannel();
		logger.debug("multiplexing to channel " + channelName);
		if (channels.get(channelName) != null) {
			channels.get(channelName).onChipsterMessage(msg);
		}
	}
}
