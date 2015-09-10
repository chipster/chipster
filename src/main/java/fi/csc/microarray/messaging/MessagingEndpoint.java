/*
 * Created on Jan 20, 2005
 *
 */
package fi.csc.microarray.messaging;

import javax.jms.Destination;
import javax.jms.JMSException;

import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.message.ChipsterMessage;

public interface MessagingEndpoint {

	/**
	 * Creates and returns a new MessagingTopic. Topics are needed for receiving and sending messages. 
	 *  
	 * @param topicName name of the topic (from a predefined enumeration)
	 */
	public abstract MessagingTopic createTopic(Topics.Name topicName,
			AccessMode accessMode) throws JMSException;

	/**
	 * Not multithread safe.
	 */
	public abstract void replyToMessage(ChipsterMessage original,
			ChipsterMessage reply) throws JMSException;

	/**
	 * Not multithread safe.
	 */
	public abstract void replyToMessage(ChipsterMessage original,
			ChipsterMessage reply, String replyChannel) throws JMSException;


	/**
	 * Closes endpoint and frees resources.
	 */
	public abstract void close() throws JMSException;

	public abstract AuthenticationRequestListener getAuthenticationListener();

	public abstract void setAuthenticationListener(
			AuthenticationRequestListener authenticationListener);

	public abstract String getSessionID();

	public abstract void setSessionID(String sessionID);

	public abstract void sendMessageToClientReplyChannel(Destination destination, ChipsterMessage msg) throws JMSException;

	/**
	 * For testing only.
	 * 
	 * @return
	 */

}