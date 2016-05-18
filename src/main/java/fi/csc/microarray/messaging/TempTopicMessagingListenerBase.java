package fi.csc.microarray.messaging;

import javax.jms.JMSException;


public abstract class TempTopicMessagingListenerBase implements TempTopicMessagingListener {

	MessagingTopic tempTopic;

	public void setTempTopic(MessagingTopic tempTopic) {
		this.tempTopic = tempTopic;
	}
	
	/**
	 * Try to delete the temp topic.
	 * 
	 */
	public void cleanUp()  {
		try {
			if (this.tempTopic != null) {
				this.tempTopic.delete();
			}
		} catch (JMSException e) {
			// nothing to do
		}
	}

}
