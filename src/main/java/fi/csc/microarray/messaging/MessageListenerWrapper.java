package fi.csc.microarray.messaging;

import javax.jms.JMSException;
import javax.jms.MapMessage;
import javax.jms.Message;
import javax.jms.MessageListener;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.util.Exceptions;

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
		
		String msgClass = "";
		MapMessage mapMessage = null;
		try {
			msgClass = msg.getStringProperty(ChipsterMessage.KEY_CLASS);
			logger.debug("message received, class is " + msgClass);
			mapMessage = (MapMessage)msg;
			ChipsterMessage chipsterMessage = (ChipsterMessage)Class.forName(msgClass).newInstance();
			chipsterMessage.unmarshal(mapMessage);
			actualListener.onChipsterMessage(chipsterMessage);
			
		} catch (Exception e) {
			logger.error("Exception when handling a message.", e);
			logger.error("message class was: " + msgClass);
			if (mapMessage != null) {
				try {
					logger.error("message first property:");
					if (mapMessage.getPropertyNames().hasMoreElements()) {
						String name = (String) mapMessage.getPropertyNames().nextElement();
						logger.error(name + " : " + mapMessage.getStringProperty(name));
					}

					logger.error("message first map object:");
					if (mapMessage.getMapNames().hasMoreElements()) {
						String name = (String) mapMessage.getMapNames().nextElement();
						logger.error(name + " : " + mapMessage.getString(name));
					}

					logger.error("message JMSType: " + mapMessage.getJMSType());
					logger.error("message JMSMessageID: " + mapMessage.getJMSMessageID());
					logger.error("message JMSDestination: " + mapMessage.getJMSDestination());
				
				
				} catch (JMSException e1) {
					logger.error(Exceptions.getStackTrace(e1));
				}
			}
		} 
	}

}
