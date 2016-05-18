package fi.csc.microarray.messaging;

public interface AuthMessagingListener extends MessagingListener {

	public void addPendingReplyListener(TempTopicMessagingListener replyListener);
}
