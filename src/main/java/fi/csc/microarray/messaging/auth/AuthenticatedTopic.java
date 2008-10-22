package fi.csc.microarray.messaging.auth;

import javax.jms.JMSException;
import javax.jms.Session;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener.Credentials;
import fi.csc.microarray.messaging.message.AuthenticationMessage;
import fi.csc.microarray.messaging.message.NamiMessage;
import fi.csc.microarray.messaging.message.AuthenticationMessage.AuthenticationOperation;


/**
 * Extended version of MessagingTopic that allows message receivers to request
 * authentication. If receivers requests authentication (no session, session expired etc.), 
 * callback is used to poll sender for credentials. 
 *
 */
public class AuthenticatedTopic extends MessagingTopic {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger
			.getLogger(AuthenticatedTopic.class);
	
	private AuthenticationRequestListener listener;
	private String sessionID = null;
	
	private MessagingListener authTopicListener = new MessagingListener() {
		public void onNamiMessage(NamiMessage msg) {
			
			try {
				if (msg instanceof AuthenticationMessage) {
					AuthenticationMessage authMsg = (AuthenticationMessage)msg;
					
					if (authMsg.isRequestForAuthentication()) {
						logger.debug("got request for authentication related to topic " + getName());

						sessionID = msg.getSessionID(); // record session for authenticating following messages

						if (listener != null) {
							Credentials credentials = listener.authenticationRequest();

							// send reply
							logger.debug("got authentication request, will send a reply, using session " + msg.getSessionID());
							AuthenticationMessage replyMsg = new AuthenticationMessage(AuthenticationOperation.LOGIN);
							replyMsg.setUsername(credentials.username);
							replyMsg.setPassword(credentials.password);
							replyMsg.setSessionID(msg.getSessionID());
							replyMsg.setReplyTo(authMsg.getReplyTo());

							sendMessage(replyMsg);						
						}
						
					} else if (authMsg.isLoginAck()) {						
						listener.authenticationSucceeded();
						
					} else {
						throw new RuntimeException("unknown authentication operation type: " + authMsg.getCommand());
					}

				} else {
					throw new RuntimeException("illegal message type received: " + msg.getClass().getSimpleName());
				}
				
			} catch (JMSException e) {
				throw new RuntimeException(e);
			}
			
		}
	};

	public AuthenticatedTopic(Session session, String topicName, Type type, AccessMode accessMode, AuthenticationRequestListener listener, MessagingEndpoint endpoint) throws JMSException {
		super(session, topicName, type, accessMode, endpoint);
		this.listener = listener;
	}
	
	/**
	 * @see MessagingTopic#sendMessage(NamiMessage)
	 */
	@Override
	public void sendMessage(NamiMessage message) throws JMSException {
		attachSessionID(message);
		super.sendMessage(message);
	}
	
	public void sendReplyableMessage(NamiMessage message, MessagingListener replyListener) throws JMSException {
		attachSessionID(message);
		logger.debug("added authentication listener to message");
		super.sendReplyableMessage(message, replyListener, authTopicListener);
	}
	
	private void attachSessionID(NamiMessage msg) {
		if (sessionID != null) {
			msg.setSessionID(sessionID);
		}
	}
}
