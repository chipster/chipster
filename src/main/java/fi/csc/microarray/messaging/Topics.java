package fi.csc.microarray.messaging;

/**
 * Helper class for handling messaging topics.
 *  
 * @author Aleksi Kallio
 *
 */
public interface Topics {
	
	/**
	 * Contains ALL topic names used in this system.
	 * 
	 * | *topic name*| *senders* | *receivers* | *description* |
	 * | test-topic | unit tests | frontend | test messaging |
	 * | request-topic | client | frontend | requests from client |
	 * | authorised-request-topic | frontend | analyser | analysis requests authorised by front end |
	 * | admin-topic | client | all | administrative requests |
	 * 
	 */
	enum Name {
		TEST_TOPIC("test-topic"),
		REQUEST_TOPIC("request-topic"),
		AUTHORISED_REQUEST_TOPIC("authorised-request-topic"),
		ADMIN_TOPIC("admin-topic"),
		MANAGER_TOPIC("job-log-topic"),
		AUTH_LOG_TOPIC("auth-log-topic")
		;

		private String topicName;

		private Name(String topicName) {
			this.topicName = topicName;
		}
		
		/**
		 * Return the actual JMS topic name.
		 */
		@Override
		public String toString() {
			return topicName;
		}
		
		
	}
	
	enum MultiplexName {
		REPLY_TO("reply-to"),
		AUTHORISE_TO("authorise-to");

		private String multiplexName;

		private MultiplexName(String multiplexName) {
			this.multiplexName = multiplexName;
		}
		
		/**
		 * Returns the actual multiplexing name.
		 */
		@Override
		public String toString() {
			return multiplexName;
		}
	}
}
