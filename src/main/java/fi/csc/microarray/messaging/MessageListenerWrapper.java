package fi.csc.microarray.messaging;

import javax.jms.MapMessage;
import javax.jms.Message;
import javax.jms.MessageListener;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.message.ChipsterMessage;

/**
 * For converting JMS-type messages into Nami type messages.
 * @author akallio
 *
 */
public class MessageListenerWrapper implements MessageListener {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger
			.getLogger(MessageListenerWrapper.class);

	private MessagingListener actualListener;
	
	public MessageListenerWrapper(MessagingListener listener) {
		actualListener = listener;
	}

	public void onMessage(Message msg) {
		
		try {
			String msgClass = msg.getStringProperty(ChipsterMessage.KEY_CLASS);
			logger.debug("message received, class is " + msgClass);
			MapMessage mapMessage = (MapMessage)msg;
			ChipsterMessage chipsterMessage = (ChipsterMessage)Class.forName(msgClass).newInstance();
			chipsterMessage.unmarshal(mapMessage);
			actualListener.onChipsterMessage(chipsterMessage);
			
		} catch (Exception e) {
			logger.error(e, e);
		} 
	}

}
