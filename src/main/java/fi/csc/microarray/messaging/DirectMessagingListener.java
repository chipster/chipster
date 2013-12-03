package fi.csc.microarray.messaging;

import fi.csc.microarray.messaging.message.ChipsterMessage;

/**
 * Like MessagingListener, but passes replyEndpoint along with the message so that receiver can simultaneusly communicate
 * through multiple endpoints.
 */
public interface DirectMessagingListener extends MessagingListener {

	/**
	 * Called when a new message is received. 
	 * 
	 * @param msg unmarshalled message
	 */
	public void onChipsterMessage(ChipsterMessage msg, MessagingEndpoint replyEndpoint);
	
	// TODO add method for stopping listening

}
