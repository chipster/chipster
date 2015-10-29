package fi.csc.microarray.messaging;

/**
 * Helper class for handling messaging topics.
 *  
 * @author Aleksi Kallio
 *
 */
public interface Topics {
	
	enum Name {
		TEST_TOPIC("test-topic"),
		REQUEST_TOPIC("request-topic"),
		AUTHORISED_REQUEST_TOPIC("authorised-request-topic"),
		FILEBROKER_TOPIC("filebroker-topic"),
		AUTHORISED_FILEBROKER_TOPIC("authorised-filebroker-topic"),
		FILEBROKER_ADMIN_TOPIC("filebroker-admin-topic"),
		COMP_ADMIN_TOPIC("comp-admin-topic"),
		JOBMANAGER_ADMIN_TOPIC("jobmanager-admin-topic"),
		JOBMANAGER_TOPIC("jobmanager-topic"),
		ADMIN_TOPIC("admin-topic"),
		JOB_LOG_TOPIC("job-log-topic"),
		FEEDBACK_TOPIC("feedback-topic"),
		AUTHORISED_FEEDBACK_TOPIC("authorised-feedback-topic"),
		AUTH_LOG_TOPIC("auth-log-topic"),
		AUTHORIZED_MANAGED_REQUEST_TOPIC("authorized-managed-request-topic")
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
