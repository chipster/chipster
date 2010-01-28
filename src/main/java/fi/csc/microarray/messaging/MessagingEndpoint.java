/*
 * Created on Jan 20, 2005
 *
 */
package fi.csc.microarray.messaging;

import java.io.InputStream;
import java.io.OutputStream;

import javax.jms.Connection;
import javax.jms.Destination;
import javax.jms.JMSException;
import javax.jms.MapMessage;
import javax.jms.Session;
import javax.jms.Topic;
import javax.jms.TopicConnectionFactory;

import org.apache.activemq.ActiveMQConnection;
import org.apache.activemq.ActiveMQConnectionFactory;
import org.apache.log4j.Logger;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.MessagingTopic.Type;
import fi.csc.microarray.messaging.auth.AuthenticatedTopic;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.NamiMessage;
import fi.csc.microarray.util.KeyAndTrustManager;
import fi.csc.microarray.util.UrlTransferUtil;


/**
 * A gateway to Chipster messaging fabric. Implements also admin features, for 
 * managing the message fabric. 
 * 
 * @author Aleksi Kallio
 *
 */
public class MessagingEndpoint implements MessagingListener {
	
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(MessagingEndpoint.class);

	/**
	 * ActiveMQ keyword for redialling.
	 */ 	
	private static final String RELIABLE_CONNECTION_SPECIFIER = "failover:"; // since ActiveMQ 4.xx, "reliable" is "failover"

	/**
	 * The broker to connect to.
	 */
	private final String brokerUrl;
	
	/**
	 * Is redialling enabled?
	 */
	private final boolean useReliable;

	private final String DEFAULT_REPLY_CHANNEL = Topics.MultiplexName.REPLY_TO.toString();
	
	private ActiveMQConnection connection;
	private MessagingTopic adminTopic = null;
	private Node master;
	private AuthenticationRequestListener authenticationListener;
	private String sessionID = null;

	/**
	 *  Created endpoint that is used as a gateway to communication fabric.
	 *  
	 * @param master callback interface to component using this endpoint
	 * @throws MicroarrayException
	 * @throws TimeoutExpiredException when broker can not be reached
	 */
	public MessagingEndpoint(Node master) throws MicroarrayException {
		this(master, null);
	}
	
	/**
	 * Endpoint constructor for components that can process authentication requests (ie. clients that can
	 * poll user for credentials).
	 * 
	 * @see #MessagingEndpoint(Node)
	 */
	public MessagingEndpoint(Node master, AuthenticationRequestListener authenticationListener) throws MicroarrayException {
		this.master = master;
		this.authenticationListener = authenticationListener;

		// configure everything
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();

		// fix URL file transfer proxy policy
		if (configuration.getBoolean("messaging", "disable-proxy")) {
			UrlTransferUtil.disableProxies();
		}

		String protocol = configuration.getString("messaging", "broker-protocol");
		String host = configuration.getString("messaging", "broker-host");
		int port = configuration.getInt("messaging", "broker-port");
		
		// check that all configs were properly filled
		if (protocol.trim().equals("") || host.trim().equals("")) {
			throw new RuntimeException("configuration error: protocol or host empty");
		}
		
		// set broker address
		useReliable = configuration.getBoolean("messaging", "use-reliable");
		brokerUrl =  protocol + "://" + host + ":" + port;
		
		// setup keystore if needed
		if ("ssl".equals(configuration.getString("messaging", "broker-protocol"))) {
			try {
				KeyAndTrustManager.initialise(
						configuration.getString("security", "keystore"),
						configuration.getString("security", "keypass").toCharArray(), 
						configuration.getString("security", "keyalias"), 
						configuration.getString("security", "master-keystore"));
			} catch (Exception e) {
				throw new MicroarrayException("could not access SSL keystore", e);
			}
		}
		
		String username;
		String password;
		try {
			username = configuration.getString("security", "username");
			if (username == null || username.trim().length() == 0) {
				throw new IllegalArgumentException("Username was not available from configuration");
			}

			password = configuration.getString("security", "password");
			if (password == null || password.trim().length() == 0) {
				throw new IllegalArgumentException("Password was not available from configuration");
			}
		} catch (Exception e) {
			throw new RuntimeException("reading authentication information failed: " + e.getMessage());
		}
		
		
		try {
			logger.info("connecting to " + brokerUrl);
			String completeBrokerUrl = brokerUrl;
			if (useReliable) {
				// tests connecting with unreliable, so that if broker is not available, 
				// we won't initiate retry sequence
				logger.debug("testing connecting to " + completeBrokerUrl);
				TopicConnectionFactory connectionFactory = new ActiveMQConnectionFactory(username, password, completeBrokerUrl);
				Connection tempConnection = connectionFactory.createTopicConnection();
				tempConnection.start();
				tempConnection.stop();
				tempConnection.close(); // it worked, we have a network connection
				
				// switch to reliable
				completeBrokerUrl = RELIABLE_CONNECTION_SPECIFIER + completeBrokerUrl;
			}
			
			// create actual reliable connection
			TopicConnectionFactory reliableConnectionFactory = new ActiveMQConnectionFactory(username, password, completeBrokerUrl);		
			connection = (ActiveMQConnection)reliableConnectionFactory.createTopicConnection();
			connection.setExceptionListener(master);
			connection.start();
			logger.debug("connected to " + completeBrokerUrl);
			
			// create admin topic
			adminTopic = createTopic(Topics.Name.ADMIN_TOPIC, AccessMode.READ_WRITE); // endpoint reacts to requests from admin-topic
			adminTopic.setListener(this);
			logger.debug("endpoint created succesfully");
		} catch (JMSException e) {
			throw new MicroarrayException("could not connect to message broker at " + brokerUrl + " (" + e.getMessage() + ")", e);
		}
	}

	/**
	 * Creates and returns a new MessagingTopic. Topics are needed for receiving and sending messages. 
	 *  
	 * @param topicName name of the topic (from a predefined enumeration)
	 */
	public MessagingTopic createTopic(Topics.Name topicName, AccessMode accessMode) throws JMSException {
		Session session = connection.createSession(false, Session.AUTO_ACKNOWLEDGE);
		return new AuthenticatedTopic(session, topicName.toString(), Type.NORMAL, accessMode, authenticationListener, this);		
	}
	
    public void replyToMessage(NamiMessage original, NamiMessage reply) throws JMSException {
    	replyToMessage(original, reply, DEFAULT_REPLY_CHANNEL);
    }

    public void replyToMessage(NamiMessage original, NamiMessage reply, String replyChannel) throws JMSException {
    	reply.setMultiplexChannel(replyChannel);
    	Destination replyToDest = original.getReplyTo();
    	replyToMessage(replyToDest, reply);
    }

    private void replyToMessage(Destination replyToDest, NamiMessage reply) throws JMSException {
		Session session = connection.createSession(false, Session.AUTO_ACKNOWLEDGE);
    	MapMessage msg = session.createMapMessage();
    	reply.marshal(msg);
    	session.createProducer(replyToDest).send(msg);
    }
    
	/**
	 * Reacts to administrative messages.
	 */
	public void onNamiMessage(NamiMessage msg) {
		try {
			CommandMessage txtMsg = (CommandMessage)msg;
			logger.debug("got admin request " + txtMsg.getCommand());
			
			if (txtMsg.getCommand().equals("ping")) {
				CommandMessage reply = new CommandMessage("ping-reply");
				reply.addParameter(master.getName());
				reply.addParameter(master.getHost());
				adminTopic.sendMessage(reply);
				logger.debug("sent ping-reply from " + master.getHost() + "/" + master.getName());
				
			} else if (txtMsg.getCommand().equals("request-load-info")) {
				if (master instanceof MonitoredNode) {
					CommandMessage reply = new CommandMessage("request-load-info-reply");
					MonitoredNode mMaster = (MonitoredNode)master;
					reply.addParameter(master.getName());
					reply.addParameter(Long.toString(mMaster.countRequestsInProcessing()));
					reply.addParameter(Long.toString(mMaster.getLastProcessingTime()));
					adminTopic.sendMessage(reply);
					logger.debug("sent request-load-info-reply from " + master.getHost() + "/" + master.getName());
				}
			}

		} catch (JMSException e) {
			logger.error(e);
		}
	}

	/**
	 * Closes endpoint and frees resources.
	 */
    public void close() throws JMSException {
    	connection.stop();
    	connection.close();	
    }
    
	public AuthenticationRequestListener getAuthenticationListener() {
		return authenticationListener;
	}
	
	public void setAuthenticationListener(AuthenticationRequestListener authenticationListener) {
		this.authenticationListener = authenticationListener;
	}

	public OutputStream createOutputStream(Topic topic) throws JMSException {
		return connection.createOutputStream(topic);
	}

	public InputStream createInputStream(Topic topic) throws JMSException {
		return connection.createInputStream(topic);
	}

	public String getSessionID() {
		return sessionID;
	}

	public void setSessionID(String sessionID) {
		this.sessionID = sessionID;
	}
}
