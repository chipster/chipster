package fi.csc.microarray.messaging;

import fi.csc.microarray.messaging.message.NamiMessage;

/**
 * Messaging fabric counterpart of JMS message listener.
 * 
 */
public interface MessagingListener {

	/**
	 * Called when a new message is received. 
	 * 
	 * @param msg unmarshalled message
	 */
	public void onNamiMessage(NamiMessage msg);
	
	// TODO add method for stopping listening

}
