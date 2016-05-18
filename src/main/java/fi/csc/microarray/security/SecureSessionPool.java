package fi.csc.microarray.security;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import org.apache.log4j.Logger;


import fi.csc.microarray.config.DirectoryLayout;

/**
 * Implementation of a session pool. Session are given cryptographically strong
 * pseudo random identifier. SecureSessionPool objects are thread-safe.
 *   
 * @author Aleksi Kallio
 *
 */
public class SecureSessionPool {
	
	private static final Logger logger = Logger.getLogger(SecureSessionPool.class);
	
	private long activeTimeout;
	private long totalTimeout;
	
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
	
	public SecureSessionPool() {
		// timeouts, hours from config file
		int activeTimeoutMinutes = DirectoryLayout.getInstance().getConfiguration().getInt("auth", "sessionRefreshTimeout");
		int totalTimeoutMinutes = DirectoryLayout.getInstance().getConfiguration().getInt("auth", "sessionMaxLifetime");

		this.activeTimeout = 1000 * 60 * activeTimeoutMinutes;
		this.totalTimeout = 1000 * 60 * totalTimeoutMinutes;
		
		logger.info("session lifetimes: active " + activeTimeoutMinutes +  "m, maximum " + totalTimeoutMinutes + "m");
	}
	
	public synchronized Session createSession() {
		String id = CryptoKey.generateRandom();
		Session session = new Session(id);
		sessions.put(id, session);
		return session;
	}
	
	public synchronized void addSession(Session session) {
		sessions.put(session.getID(), session);
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
		return (System.currentTimeMillis() - session.getLastUseTimestamp()) > activeTimeout ||
			(System.currentTimeMillis() - session.getCreationTimestamp()) > totalTimeout;
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
