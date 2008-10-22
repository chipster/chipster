package fi.csc.microarray.messaging;

import javax.jms.MapMessage;
import javax.jms.Message;
import javax.jms.MessageListener;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.message.NamiMessage;

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
			String msgClass = msg.getStringProperty(NamiMessage.KEY_CLASS);
			logger.debug("message received, class is " + msgClass);
			MapMessage mapMessage = (MapMessage)msg;
			NamiMessage namiMessage = (NamiMessage)Class.forName(msgClass).newInstance();
			namiMessage.unmarshal(mapMessage);
			actualListener.onNamiMessage(namiMessage);
			
		} catch (Exception e) {
			logger.error(e, e);
		} 
	}

}
