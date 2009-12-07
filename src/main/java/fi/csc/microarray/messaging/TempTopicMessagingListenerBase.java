package fi.csc.microarray.messaging;

import javax.jms.JMSException;

import org.apache.log4j.Logger;


public abstract class TempTopicMessagingListenerBase implements TempTopicMessagingListener {

	private static final Logger logger = Logger.getLogger(TempTopicMessagingListenerBase.class);
	
	MessagingTopic tempTopic;

	public void setTempTopic(MessagingTopic tempTopic) {
		this.tempTopic = tempTopic;
	}
	
	public void cleanUp()  {
		try {
			this.tempTopic.delete();
		} catch (JMSException e) {
			logger.warn("Could not delete temporary topic", e);
		}
	}
}
