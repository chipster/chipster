/*
 * Created on Feb 24, 2005
 *
*/
package fi.csc.microarray.messaging.message;

import java.util.UUID;

import javax.jms.Destination;
import javax.jms.JMSException;
import javax.jms.MapMessage;

import org.apache.log4j.Logger;

/**
 * The base class for all messages handled by the NAMI fabric.
 * 
 * @author hupponen, akallio
 *
 */
public abstract class ChipsterMessage {

	private static final Logger logger = Logger.getLogger(ChipsterMessage.class);
	
	public static final String KEY_MESSAGE_ID = "message-id";
	public static final String KEY_CLASS= "class";	
	public static final String KEY_USERNAME = "username";
	public static final String KEY_SESSION_ID = "session-id";
	public static final String KEY_MULTIPLEX_CHANNEL = "multiplex-channel";
    
	private String messageID = UUID.randomUUID().toString();
	private Destination replyTo = null;
	private String username;
	private String sessionID;
	private String multiplexChannel;
	private String jmsMessageID;
	
	/**
	 * Converts JMS MapMessage into ChipsterMessage. Class extenders, see marshal-method.
	 * @see #marshal(MapMessage)
	 */
	public void unmarshal(MapMessage from) throws JMSException {
		this.messageID = from.getStringProperty(KEY_MESSAGE_ID);
		this.replyTo = from.getJMSReplyTo();
		this.username = from.getStringProperty(KEY_USERNAME);
		this.multiplexChannel = from.getStringProperty(KEY_MULTIPLEX_CHANNEL);
		this.sessionID = from.getStringProperty(KEY_SESSION_ID);
		this.jmsMessageID = from.getJMSMessageID();
	}
	
	
	/**
	 * <p>Extenders should first call superclass' marshal-method. 
	 * Marshalling should use JMS properties (e.g. setStringProperty)
	 * for non-payload data and JMS MapMessage values (e.g. setString)
	 * for payload data ie. small data as properties and big data
	 * as values. There are no hard reasons for this, but this way we 
	 * give hints to JMS implementation on how to handle our data and
	 * also make debugging easier, as properties are printed out to logs, 
	 * but values are not. For this reason, sensitive data should be stored
	 * as values. Also only properties can be accessed with message selectors.</p>
	 *  
	 * @param to
	 * @throws JMSException
	 */
	public void marshal(MapMessage to) throws JMSException {
		to.setStringProperty(KEY_MESSAGE_ID, messageID);
		to.setStringProperty(KEY_CLASS, this.getClass().getCanonicalName());
		to.setStringProperty(KEY_USERNAME, username);
		to.setJMSReplyTo(replyTo);
		to.setStringProperty(KEY_MULTIPLEX_CHANNEL, multiplexChannel);
		to.setStringProperty(KEY_SESSION_ID, sessionID);
	}
	
	/**
	 * Returns the unique identifier of the message. Set by front-end or other 
	 * party capable of guarenteeing uniqueness.
	 */
	public String getMessageID() {
		return messageID;
	}
		
	/**
	 * Returns the Destination used for replying to this message.
	 */
	public Destination getReplyTo() {
		return replyTo;
	}

	/**
	 * Reply-to is set by the messaging system. Use MessagingTopic.sendReplyableMessage
	 * for sending a replyable message.
	 *  
	 * @see #getReplyTo()
	 * @see fi.csc.microarray.messaging.MessagingTopic#sendReplyableMessage(ChipsterMessage, MessagingListener, MessagingListener)
	 */
	public void setReplyTo(Destination replyTo) {
		this.replyTo = replyTo;
	}
	
	public String toString() {
		return "MESSAGE message-id: " + messageID + ", class: "+ this.getClass().getSimpleName();  
	}

	/**
	 * Returns the username part of credentials associated with this message.
	 */
	public String getUsername() {
		return username;
	}

	/**
	 * @see #getUsername()
	 */
	public void setUsername(String username) {
		this.username = username;
	}

	public String getSessionID() {
		return sessionID;
	}

	public void setSessionID(String sessionID) {
		this.sessionID = sessionID;
	}

	/**
	 * Multiplex channel is used to transfer multiple topics across one physical topic. 
	 */
	public String getMultiplexChannel() {
		return multiplexChannel;
	}

	/**
	 * @see #getMultiplexChannel()
	 */
	public void setMultiplexChannel(String multiplexChannel) {
		this.multiplexChannel = multiplexChannel;
	}


	public String getJmsMessageID() {
		return jmsMessageID;
	}

	public void handleException(Exception e) throws JMSException {
		logger.error(e);
		throw new JMSException(e.getMessage()); // converting URL related errors to JMS errors is kind of strange...
	}


}
