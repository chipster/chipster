package fi.csc.microarray.messaging.auth;

import java.util.LinkedList;
import java.util.List;

import javax.jms.JMSException;
import javax.jms.Session;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.AuthMessagingListener;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.TempTopicMessagingListener;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener.Credentials;
import fi.csc.microarray.messaging.message.AuthenticationMessage;
import fi.csc.microarray.messaging.message.AuthenticationMessage.AuthenticationOperation;
import fi.csc.microarray.messaging.message.ChipsterMessage;


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
	
	private AuthMessagingListener authMessagingListener = new AuthMessagingListener() {
		private List<TempTopicMessagingListener> pendingReplyListeners = new LinkedList<TempTopicMessagingListener>();
		
		public void onChipsterMessage(ChipsterMessage msg) {
			
			try {
				if (msg instanceof AuthenticationMessage) {
					AuthenticationMessage authMsg = (AuthenticationMessage)msg;
					
					if (authMsg.isRequestForAuthentication()) {
						logger.debug("got request for authentication related to topic " + getName());

				
						if (listener != null) {
							Credentials credentials = listener.authenticationRequest();
							
							// cancel from dialog, signal pending reply listeners 
							if (credentials == null) {
								for (TempTopicMessagingListener pendingReplyListener : pendingReplyListeners) {
									pendingReplyListener.cancel();
								}
								return;
							}
							
							// send reply
							logger.debug("got authentication request, will send a reply, using session " + msg.getSessionID());
							AuthenticationMessage replyMsg = new AuthenticationMessage(AuthenticationOperation.LOGIN);
							replyMsg.setUsername(credentials.username);
							replyMsg.setPassword(credentials.password);
							replyMsg.setSessionID(msg.getSessionID());
							replyMsg.setReplyTo(authMsg.getReplyTo());
							
							// FIXME sometimes reply to this send (ack) is not received as reply temp topic
							// gets deleted before that
							// Should send replyable message.
							sendMessage(replyMsg);
						}
						
					} else if (authMsg.isLoginAck()) {						
						fi.csc.microarray.client.Session.getSession().setUsername(authMsg.getUsername());
						getEndpoint().setSessionID(authMsg.getSessionID()); // record session for authenticating following messages
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

		@Override
		public void addPendingReplyListener(TempTopicMessagingListener pendingReplyListener) {
			this.pendingReplyListeners.add(pendingReplyListener);
		}
	};

	public AuthenticatedTopic(Session session, String topicName, Type type, AccessMode accessMode, AuthenticationRequestListener listener, MessagingEndpoint endpoint) throws JMSException {
		super(session, topicName, type, accessMode, endpoint);
		this.listener = listener;
	}
	
	/**
	 * @see MessagingTopic#sendMessage(ChipsterMessage)
	 */
	@Override
	public void sendMessage(ChipsterMessage message) throws JMSException {
		// attach session id to messages other than login messages
		if (!(message instanceof AuthenticationMessage && ((AuthenticationMessage)message).isLogin())) {
			attachSessionID(message);
		}
		super.sendMessage(message);
	}
	
	public void sendReplyableMessage(ChipsterMessage message, TempTopicMessagingListener replyListener) throws JMSException {
		// attach session id to messages other than login messages
		if (!(message instanceof AuthenticationMessage && ((AuthenticationMessage)message).isLogin())) {
			attachSessionID(message);
		}
		logger.debug("added authentication listener to message");
		super.sendReplyableMessage(message, replyListener, authMessagingListener);
	}
	
	private void attachSessionID(ChipsterMessage msg) {
		String sessionID = getEndpoint().getSessionID();
		if (sessionID != null) {
			msg.setSessionID(sessionID);
		}
	}
}
