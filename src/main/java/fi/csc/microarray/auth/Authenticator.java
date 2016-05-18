package fi.csc.microarray.auth;

import java.util.Arrays;

import javax.jms.Destination;
import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.messaging.JMSMessagingEndpoint;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.message.AuthenticationMessage;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.AuthenticationMessage.AuthenticationOperation;
import fi.csc.microarray.security.SecureSessionPool;
import fi.csc.microarray.security.SecureSessionPool.Session;
import fi.csc.microarray.service.KeepAliveShutdownHandler;
import fi.csc.microarray.service.ShutdownCallback;
import fi.csc.microarray.util.SystemMonitorUtil;

/**
 * @author Aleksi Kallio
 *
 */
public class Authenticator extends NodeBase implements ShutdownCallback {
    
	/**
	 * Logger for this class
	 */
	private static Logger logger = null;

	/**
	 * Logger for security issues
	 */
	private static Logger securityLogger = null;

	/**
	 *  Logger for message logging
	 */
	private static Logger messageLogger = null;
	
	private SecureSessionPool pendingSessions;
	private SecureSessionPool validSessions;
	private MessagingEndpoint endpoint;
	private MessagingTopic authorisedTopic;
	private MessagingTopic authorisedFilebrokerTopic;
	private MessagingTopic authorisedFeedbackTopic;

	private MessagingTopic testTopic;

	private TestListener testListener;
	private AuthenticationProvider authenticationProvider; 
	
	public Authenticator(String configURL) throws Exception {
		
		// initialise dir and logging
		DirectoryLayout.initialiseServerLayout(Arrays.asList(new String[] {"auth"}), configURL);
		logger = Logger.getLogger(Authenticator.class);
		securityLogger = Logger.getLogger("security.frontend");
		messageLogger = Logger.getLogger("messages.frontend");
		
		// initialise session pools
		// FIXME set proper lifetimes for pending sessions
		pendingSessions = new SecureSessionPool();
		validSessions = new SecureSessionPool();
		
		// initialise communications
		this.endpoint = new JMSMessagingEndpoint(this);

		// create authorised topics
		authorisedTopic = endpoint.createTopic(Topics.Name.AUTHORISED_REQUEST_TOPIC, AccessMode.WRITE);       
		authorisedFilebrokerTopic = endpoint.createTopic(Topics.Name.AUTHORISED_FILEBROKER_TOPIC, AccessMode.WRITE);       
		authorisedFeedbackTopic = endpoint.createTopic(Topics.Name.AUTHORISED_FEEDBACK_TOPIC, AccessMode.WRITE);       

		
		// create non-authorised topics
		RequestListener jobListener = new RequestListener(authorisedTopic);
		MessagingTopic requestTopic = endpoint.createTopic(Topics.Name.REQUEST_TOPIC, AccessMode.READ);		
		requestTopic.setListener(jobListener);		
		RequestListener filebrokerListener = new RequestListener(authorisedFilebrokerTopic);
		MessagingTopic filebrokerTopic = endpoint.createTopic(Topics.Name.FILEBROKER_TOPIC, AccessMode.READ);		
		filebrokerTopic.setListener(filebrokerListener);
		RequestListener feedbackListener = new RequestListener(authorisedFeedbackTopic);
		MessagingTopic feedbackTopic = endpoint.createTopic(Topics.Name.FEEDBACK_TOPIC, AccessMode.READ);		
		feedbackTopic.setListener(feedbackListener);

		
		// create test-topic
		testTopic = endpoint.createTopic(Topics.Name.TEST_TOPIC, AccessMode.READ_WRITE);		
		this.testListener = new TestListener();
		testTopic.setListener(testListener);
		
		// initialise JAAS authentication
		authenticationProvider = new JaasAuthenticationProvider();
		
		// create keep-alive thread and register shutdown hook
		KeepAliveShutdownHandler.init(this);
		
		logger.info("authenticator is up and running [" + ApplicationConstants.VERSION + "]");
		logger.info("[mem: " + SystemMonitorUtil.getMemInfo() + "]");
	}
	
	private class RequestListener implements MessagingListener {

		private static final String KEY_PENDING_MESSAGE = "pending-message";
		private static final String KEY_USERNAME = "username";
		
		private MessagingTopic routeTo;

		public RequestListener(MessagingTopic routeTo) {
			this.routeTo = routeTo;
		}

		/**
		 * Two step processing: authenticate, then route.
		 */
		public void onChipsterMessage(ChipsterMessage msg) {
			try {
				logger.debug("starting to process " + msg);

				// pick login messages
				if (msg instanceof AuthenticationMessage && ((AuthenticationMessage)msg).isLogin()) {
					login((AuthenticationMessage)msg);
					return;
				}
				
				// try to load existing session
				Session session = null;
				if (msg.getSessionID() != null) {
					String id = msg.getSessionID();
					if (validSessions.getSession(id) != null) {
						session = validSessions.getSession(id);
						logger.debug("message " + msg.getMessageID() + " had a proper session " + session.getID());
					}					
				}

				// message had no valid session --> request authentication
				if (session == null) {
					logger.debug("message " + msg + " has no session, requires authentication");
					// 1.2.a. no existing session (possibly due to expiration)

					// create session and request authentication
					session = pendingSessions.createSession();
					logger.debug("created new session, pending sessions pool size is now " + pendingSessions.size());
					session.putParameter(KEY_PENDING_MESSAGE, msg);
					
					requestAuthentication(msg, session.getID().toString());
					return;
				} 
				
				// message had a valid session
				else {					
					
					// logout message (login has been taken care before)
					if (msg instanceof AuthenticationMessage) {
						AuthenticationMessage authMsg = (AuthenticationMessage)msg;
						if (authMsg.isLogout()) {
							logger.debug("message " + msg.getMessageID() + " is a logout");
							validSessions.removeSession(session);
							return;
						} else {
							throw new IllegalArgumentException("unknown authentication message: " + authMsg);
						}
					} 
					
					// all other messages
					else {
						routeMessage(msg, session);
						session.touch(); // session succesfully used, so update timestamp
					}
				}
				
			} catch (JMSException e) {
				logger.error(e);
			} catch (AuthorisationException e) {
				securityLogger.info("authorisation failed: " +  e.getMessage());
				logger.info(e);
			} 
		}

		
		private void login(AuthenticationMessage authMsg) throws JMSException, AuthorisationException {
			messageLogger.debug("message " + authMsg.getMessageID() + " is a login");

			// get session id from message
			String sessionId = authMsg.getSessionID();
			if (sessionId == null || sessionId.isEmpty()) {
				logger.warn("got login message with invalid session id: " + sessionId);
				return;
			}
			Session session = pendingSessions.getSession(sessionId);
			
			// get session
			if (session == null) {
				logger.warn("pending session " + sessionId + " not found");
				return;
			} 
			
			// authenticate with username/password ok
			if (authenticationProvider.authenticate(authMsg.getUsername(), authMsg.getPassword().toCharArray())) {
				
				pendingSessions.removeSession(session);
				session.putParameter(KEY_USERNAME, authMsg.getUsername());
				validSessions.addSession(session);
				securityLogger.info("authenticated user " + authMsg.getUsername() + " (auth. message JMS id was " + authMsg.getJmsMessageID() + ")");

				// ack username to client
				try {
					ackLogin(authMsg, session.getID().toString(), true);
				} catch (Exception e) {
					logger.warn("could not send acknowledge message for " + authMsg.getUsername());
				}

				// send possible pending message
				ChipsterMessage pendingMessage = (ChipsterMessage)session.getParameter(KEY_PENDING_MESSAGE);
				if (pendingMessage != null) {
					// reply to first login message
					if (pendingMessage instanceof CommandMessage && 
							CommandMessage.COMMAND_LOGIN.equals(((CommandMessage)pendingMessage).getCommand())) {
						sendReplyToFirstLogin((CommandMessage)pendingMessage);
					}
					
					// normal messages
					else {
						routeMessage(pendingMessage, session);
					}
				}
				session.touch();
			} 
			
			// authentication with username/password failed
			else {
				securityLogger.info("illegal username/password (user " + authMsg.getUsername()  + ", auth. message JMS id was " + authMsg.getJmsMessageID() + ")");
				ackLogin(authMsg, session.getID().toString(), false);
				return;
			}
		}

		private void sendReplyToFirstLogin(CommandMessage loginMessage) {
			AuthenticationMessage reply = new AuthenticationMessage(AuthenticationOperation.LOGIN_SUCCEEDED);
			logger.debug("sending reply to first login message");
			try {
				endpoint.replyToMessage(loginMessage, reply, Topics.MultiplexName.REPLY_TO.toString());
			} catch (JMSException e) {
				logger.warn("failed to send reply to first login");
			}
		}

		private void routeMessage(ChipsterMessage message, Session session) throws JMSException {
			// sanity check for username;
			String username = (String)session.getParameter(KEY_USERNAME);
			if (username == null || username.equals("")) {
				logger.warn("not routing a message with null or empty username");
				Session sessionToBeRemoved = validSessions.getSession(session.getID()); 
				if (sessionToBeRemoved != null) {
					validSessions.removeSession(sessionToBeRemoved);	
				}
				return;
			}
			
			// prepare
			message.setUsername(username);
			message.setSessionID(null);
			
			// send
			messageLogger.info(message);
			routeTo.sendMessage(message);
	
		}
		
		
		private void ackLogin(ChipsterMessage loginMessage, String sessionID, boolean succeeded) throws JMSException, AuthorisationException {
			
			AuthenticationOperation operation = succeeded ? AuthenticationMessage.AuthenticationOperation.LOGIN_SUCCEEDED : AuthenticationMessage.AuthenticationOperation.LOGIN_FAILED; 
			AuthenticationMessage ackMessage = new AuthenticationMessage(operation);
			ackMessage.setSessionID(sessionID); // needed both for success and fail
			ackMessage.setReplyTo(loginMessage.getReplyTo());
			if (succeeded) {
				ackMessage.setUsername(loginMessage.getUsername());
			}
			endpoint.replyToMessage(loginMessage, ackMessage, Topics.MultiplexName.AUTHORISE_TO.toString());
		}
		
		private void requestAuthentication(ChipsterMessage msg, String sessionID) throws JMSException, AuthorisationException {
			AuthenticationMessage request = new AuthenticationMessage(AuthenticationMessage.AuthenticationOperation.REQUEST);
			request.setSessionID(sessionID);
			request.setReplyTo(msg.getReplyTo());
			logger.debug("requesting authentication for " + msg.getMessageID() + " under the session " + sessionID);
			endpoint.replyToMessage(msg, request, Topics.MultiplexName.AUTHORISE_TO.toString());
		}
	}
		
	private class TestListener implements MessagingListener {
		public void onChipsterMessage(ChipsterMessage msg) {
			logger.debug("got message on test-topic.");
			try {

				Destination dest = msg.getReplyTo();
				if (dest != null) {
					CommandMessage replyMessage = new CommandMessage("OK");
					endpoint.replyToMessage(msg, replyMessage);
				} else {
					logger.debug("ReplyTo is null, so no reply was sent.");
				}
			} catch (Exception e) {
				logger.error(e);
			}
		}
	}

	public String getName() {
		return "authenticator";
	}

	public void shutdown() {
		logger.info("shutdown requested");

		// close messaging endpoint
		try {
			this.endpoint.close();
		} catch (JMSException e) {
			logger.error("closing messaging endpoint failed", e);
		}

		logger.info("shutting down");
	}
}
