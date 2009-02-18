package fi.csc.microarray.security;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Implementation of a session pool. Session are given cryptographically strong
 * pseudo random identifier. SecureSessionPool objects are thread-safe.
 *   
 * @author Aleksi Kallio
 *
 */
public class SecureSessionPool {
	
	private static final long ACTIVE_TIMEOUT_IN_MILLIS = 1*60*60*1000; // 1 hour
	private static final long TOTAL_TIMEOUT_IN_MILLIS = 14*60*60*1000; // 14 hours
	
	private Map<String, Session> sessions = new HashMap<String, Session>();
	
	public class Session {
				
		private long creationTimestamp;
		private long lastUseTimestamp;
		
		private String id;
		private Map<String, Object> parameters = 
			Collections.synchronizedMap(new HashMap<String, Object>());

		public Session(String id) {
			this.id = id;
			this.creationTimestamp = System.currentTimeMillis();
			touch();
		}
		
		public void putParameter(String key, Object parameter) {
			parameters.put(key, parameter);
		}
		
		public Object getParameter(String key) {
			return parameters.get(key);
		}
		
		public void removeParameter(String key) {
			parameters.remove(key);
		}
		
		public String getID() {
			return id;
		}
		
		@Override
		public String toString() {
			return parameters.toString();
		}

		/**
		 * Update last-use timestamp.
		 *
		 */
		public void touch() {
			this.lastUseTimestamp = System.currentTimeMillis();
		}

		private long getLastUseTimestamp() {
			return lastUseTimestamp;
		}
		
		private long getCreationTimestamp() {
			return creationTimestamp;
		}
	}
	
	public synchronized Session createSession() {
		String id = CryptoKey.generateRandom();
		Session session = new Session(id);
		sessions.put(id, session);
		return session;
	}
	
	public synchronized Session getSession(String id) {
		Session session = sessions.get(id);
		if (session != null && isExpired(session)) {
			// session is old
			sessions.remove(session);
			return null;
		}
		return session;
	}
	
	private boolean isExpired(Session session) {
		return (System.currentTimeMillis() - session.getLastUseTimestamp()) > ACTIVE_TIMEOUT_IN_MILLIS ||
			(System.currentTimeMillis() - session.getCreationTimestamp()) > TOTAL_TIMEOUT_IN_MILLIS;
	}
	
	public synchronized void removeSession(Session session) {
		sessions.remove(session.id);
	}

	public int size() {
		return sessions.size();
	}

	public static void main(String[] args) {
		for (int i = 0; i < 10; i++) {
			SecureSessionPool ssp = new SecureSessionPool();
			System.out.println(ssp.createSession().id);
		}
	}
}
