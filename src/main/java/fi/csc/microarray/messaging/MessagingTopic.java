/*
 * Created on Jan 20, 2005
 *
 */
package fi.csc.microarray.messaging;

import javax.jms.JMSException;
import javax.jms.MapMessage;
import javax.jms.MessageConsumer;
import javax.jms.MessageProducer;
import javax.jms.Session;
import javax.jms.TemporaryTopic;
import javax.jms.Topic;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.message.ChipsterMessage;

/**
 * One topic ("channel") in the messaging fabric. Topics are publish-subscribe
 * type of communication. Used for sending and receiving messages.
 *  
 * @author akallio
 */
public class MessagingTopic {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(MessagingTopic.class);

	public static enum Type {
		NORMAL,
		TEMPORARY,
		NAMED_TEMPORARY,
		MOCK;
	}
	
	public static enum AccessMode {
		READ,
		WRITE,
		READ_WRITE;
	}
	
	
	private MessageConsumer consumer;
	private Session session;
	private Topic topic;
	private MessagingEndpoint endpoint;
	
	public MessagingTopic(Session session, String topicName, Type type, AccessMode accessMode, MessagingEndpoint endpoint) throws JMSException {
		this.session = session;
		
		switch (type) {
		case NORMAL:
			this.topic = session.createTopic(topicName);
			break;
		case TEMPORARY:
			this.topic = session.createTemporaryTopic();
			break;
		case MOCK:
			// initialise nothing
			break;
		default:
			throw new RuntimeException("unknown type " + type);
		}
		
		// create consumer only if we need to read from the topic
		if (accessMode == AccessMode.READ || accessMode == AccessMode.READ_WRITE) {
			this.consumer = session.createConsumer(topic);	
		}
		this.endpoint = endpoint;
	}

	

	/**
	 * Not multithread safe.
	 */
	protected void sendReplyableMessage(ChipsterMessage message, TempTopicMessagingListener replyListener, AuthMessagingListener authenticationListener) throws JMSException {
		MessagingTopic tempTopic = new MessagingTopic(session, null, Type.TEMPORARY, AccessMode.READ_WRITE, endpoint);
		
		MultiplexingMessagingListener plexer = new MultiplexingMessagingListener();
		plexer.addChannel(Topics.MultiplexName.REPLY_TO.toString(), replyListener);
		if (authenticationListener != null) {
			plexer.addChannel(Topics.MultiplexName.AUTHORISE_TO.toString(), authenticationListener);
			authenticationListener.addPendingReplyListener(replyListener);
		}
		tempTopic.setListener(plexer);
		
		replyListener.setTempTopic(tempTopic);
		message.setReplyTo(tempTopic.topic);
		sendMessage(message);
		
	}
	
	/**
	 * Sends the message and creates a temporary topic for replying.
	 * Not multithread safe.
	 * 
	 * @param replyListener receives replies (if any) through hidden temporary topic
	 */
	public void sendReplyableMessage(ChipsterMessage message, TempTopicMessagingListener replyListener) throws JMSException {
		sendReplyableMessage(message, replyListener, null);
	}
	
	/**
	 * The basic message sending method. Sends a message without reply possibility.
	 * Not multithread safe.
	 */
	public void sendMessage(ChipsterMessage message) throws JMSException {

		// log
		logger.debug("sending " + message);
		
		// marshal message to MapMessage
		MapMessage mapMessage = session.createMapMessage();
		message.marshal(mapMessage);
		
		MessageProducer producer = null;
		try {
			producer = session.createProducer(topic);
			producer.send(mapMessage);
		} finally {
			try {
				producer.close();
			} catch (Exception e) {
			}
		}
	}

	/**
	 * Returns the JMS-name of the topic.
	 */
	public String getName() throws JMSException {
		return topic.getTopicName();
	}

	/**
	 * Registers a listener that receives all messages sent in this fabric
	 * to this topic. There can be only one listener. 
	 */
	public void setListener(MessagingListener listener) throws JMSException {
		if (consumer != null) {
			consumer.setMessageListener(new MessageListenerWrapper(listener));
		} else {
			throw new IllegalStateException("Topic was created as write only");
		}
	}
	
	/**
	 * Unregisters a listener.
	 */
	public void removeListener() throws JMSException {
		if (consumer != null) {
			consumer.setMessageListener(null);
		} else {
			throw new IllegalStateException("Topic was created as write only");
		}
	}

	public MessagingEndpoint getEndpoint() {
		return endpoint;
	}
	
	/**
	 * 
	 * 
	 * @throws JMSException
	 */
	public void delete() throws JMSException {
		if (this.topic instanceof TemporaryTopic) {
			this.consumer.close();
			((TemporaryTopic)topic).delete();
		}
	}
	
	public Topic getJMSTopic() {
		return this.topic;
	}
}
