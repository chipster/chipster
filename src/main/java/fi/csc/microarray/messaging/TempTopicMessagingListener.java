package fi.csc.microarray.messaging;

public interface TempTopicMessagingListener extends MessagingListener {

	public void setTempTopic(MessagingTopic tempTopic);
	
	public void cleanUp();
	
	public void cancel();
}
