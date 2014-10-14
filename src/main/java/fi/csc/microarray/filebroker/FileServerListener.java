package fi.csc.microarray.filebroker;

import java.util.List;

import fi.csc.microarray.messaging.MessagingEndpoint;

/**
 * Listen for FileServer operations.
 * 
 * @author klemela
 */
public abstract class FileServerListener {
	
	public abstract void listen(Event e);
	
	public static abstract class Event {
		private String username;
		private String uuid;
		private MessagingEndpoint endpoint;
		
		public Event(String uuid, String username, MessagingEndpoint endpoint) {
			this.uuid = uuid;
			this.username = username;
			this.endpoint = endpoint;			
		}
		
		public String getUsername() {
			return username;
		}

		public String getUuid() {
			return uuid;
		}

		public MessagingEndpoint getEndpoint() {
			return endpoint;
		}		
	}
	
	public static abstract class StoreSessionEvent extends Event {
		private String sessionName;
		private List<String> fileIds;
		
		public StoreSessionEvent(String username, String sessionName, String uuid, List<String> fileIds, MessagingEndpoint endpoint) {
			super(uuid, username, endpoint);
			this.sessionName = sessionName;
			this.fileIds = fileIds;
		}

		public String getSessionName() {
			return sessionName;
		}
		public List<String> getFileIds() {
			return fileIds;
		}
	}
		
	public static class AfterStoreSessionReply extends StoreSessionEvent {

		public AfterStoreSessionReply(String username, String sessionName, String uuid, List<String> fileIds, MessagingEndpoint endpoint) {
			super(username, sessionName, uuid, fileIds, endpoint);
		}
	}
		
	public static class BeforeRemoveSession extends Event {

		public BeforeRemoveSession(String uuid, String username, MessagingEndpoint endpoint) {
			super(uuid, username, endpoint);
		}		
	}	
}
