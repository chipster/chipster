/*
 * Created on Jan 20, 2005
 *
 */
package fi.csc.microarray.messaging;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.security.KeyStoreException;
import java.security.NoSuchAlgorithmException;
import java.security.cert.CertificateException;

import javax.jms.Connection;
import javax.jms.Destination;
import javax.jms.JMSException;
import javax.jms.MapMessage;
import javax.jms.Session;
import javax.net.ssl.SSLHandshakeException;

import org.apache.activemq.ActiveMQConnection;
import org.apache.activemq.ActiveMQConnectionFactory;
import org.apache.activemq.ActiveMQSslConnectionFactory;
import org.apache.log4j.Logger;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.MessagingTopic.Type;
import fi.csc.microarray.messaging.auth.AuthenticatedTopic;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.util.KeyAndTrustManager;
import fi.csc.microarray.util.UrlTransferUtil;


/**
 * A gateway to Chipster messaging fabric. Implements also admin features, for 
 * managing the message fabric. 
 * 
 * @author Aleksi Kallio
 *
 */
public class JMSMessagingEndpoint implements MessagingEndpoint, MessagingListener {
	
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(JMSMessagingEndpoint.class);

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
	public JMSMessagingEndpoint(Node master) throws MicroarrayException {
		this(master, null, false);
	}
	
	/**
	 * Endpoint constructor for components that can process authentication requests (ie. clients that can
	 * poll user for credentials).
	 * 
	 * @see #MessagingEndpoint(Node)
	 */
	public JMSMessagingEndpoint(Node master, AuthenticationRequestListener authenticationListener, boolean useUnreliableAtStartup) throws MicroarrayException {
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
				
				KeyAndTrustManager.initialiseTrustStore();
			} catch (Exception e) {
				throw new MicroarrayException("could not access SSL keystore", e);
			}
		}
		
		String username = getUsername(configuration);
		String password = getPassword(configuration);
				
		try {
			logger.info("connecting to " + brokerUrl);
			String completeBrokerUrl = brokerUrl;
			if (useUnreliableAtStartup) {
				// tests connecting with unreliable, so that if broker is not available, 
				// we won't initiate retry sequence
				logger.debug("testing connecting to " + completeBrokerUrl);
				ActiveMQConnectionFactory connectionFactory = createConnectionFactory(username, password, completeBrokerUrl);
				Connection tempConnection = connectionFactory.createTopicConnection();
				tempConnection.start();
				tempConnection.stop();
				
				try {
					tempConnection.close(); // it worked, we have a network connection
				} catch (Exception e) {
					logger.warn("got exception when closing test connection");
				}
			}

			// switch to reliable
			if (useReliable) {
				completeBrokerUrl = RELIABLE_CONNECTION_SPECIFIER + completeBrokerUrl;
			}
			
			// create actual reliable connection
			ActiveMQConnectionFactory reliableConnectionFactory = createConnectionFactory(username, password, completeBrokerUrl);
			connection = (ActiveMQConnection)reliableConnectionFactory.createTopicConnection();
			connection.setExceptionListener(master);
			connection.start();
			logger.debug("connected to " + completeBrokerUrl);
			
			// create admin topic
			adminTopic = createTopic(Topics.Name.ADMIN_TOPIC, AccessMode.READ_WRITE); // endpoint reacts to requests from admin-topic
			adminTopic.setListener(this);
			logger.debug("endpoint created succesfully");
			
		} catch (JMSException e) {			
			if (e.getCause() instanceof SSLHandshakeException) {								
				throw new MicroarrayException("server identity cannot be verified or other SSL error when connecting to " + brokerUrl, e);
			} else {
				throw new MicroarrayException("could not connect to message broker at " + brokerUrl, e);
			}
		}
	}
	
	private static String getPassword(Configuration configuration) {

		try {		
			String password = configuration.getString("security", "password");
			if (password == null || password.trim().length() == 0) {
				throw new IllegalArgumentException("Password was not available from configuration");
			}
			return password;
		} catch (Exception e) {
			throw new RuntimeException("reading authentication information failed: " + e.getMessage());
		}
	}

	private static String getUsername(Configuration configuration) {
		
		try {
			String username = configuration.getString("security", "username");
			if (username == null || username.trim().length() == 0) {
				throw new IllegalArgumentException("Username was not available from configuration");
			}
			return username;
		} catch (Exception e) {
			throw new RuntimeException("reading authentication information failed: " + e.getMessage());
		}
	}

	public static String getClientTruststore() throws NoSuchAlgorithmException, CertificateException, FileNotFoundException, KeyStoreException, IOException {
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();		
		return KeyAndTrustManager.getClientTrustStore(configuration, getPassword(configuration));
	}

	private ActiveMQConnectionFactory createConnectionFactory(String username, String password, String completeBrokerUrl) {

		ActiveMQSslConnectionFactory reliableConnectionFactory = new ActiveMQSslConnectionFactory();
		
		// dummy trust manager
//		reliableConnectionFactory.setKeyAndTrustManagers(null, new TrustManager[] {new X509TrustManager() {
//
//			
//			public void checkClientTrusted(X509Certificate[] chain, String authType) throws CertificateException {
//			}
//
//			public void checkServerTrusted(X509Certificate[] chain, String authType) throws CertificateException {
//			}
//
//			public X509Certificate[] getAcceptedIssuers() {
//				return null;
//			}}}, new SecureRandom());
//		}
		
		reliableConnectionFactory.setUserName(username);
		reliableConnectionFactory.setPassword(password);
		reliableConnectionFactory.setBrokerURL(completeBrokerUrl);
		

		reliableConnectionFactory.setWatchTopicAdvisories(false);
		return reliableConnectionFactory;
	}

	/* (non-Javadoc)
	 * @see fi.csc.microarray.messaging.MessagingEndpointIntrfc#createTopic(fi.csc.microarray.messaging.Topics.Name, fi.csc.microarray.messaging.MessagingTopic.AccessMode)
	 */
	@Override
	public MessagingTopic createTopic(Topics.Name topicName, AccessMode accessMode) throws JMSException {
		Session session = connection.createSession(false, Session.AUTO_ACKNOWLEDGE);
		return new AuthenticatedTopic(session, topicName.toString(), Type.NORMAL, accessMode, authenticationListener, this);		
	}
	

	/* (non-Javadoc)
	 * @see fi.csc.microarray.messaging.MessagingEndpointIntrfc#replyToMessage(fi.csc.microarray.messaging.message.ChipsterMessage, fi.csc.microarray.messaging.message.ChipsterMessage)
	 */
    @Override
	public void replyToMessage(ChipsterMessage original, ChipsterMessage reply) throws JMSException {
    	replyToMessage(original, reply, DEFAULT_REPLY_CHANNEL);
    }

	/* (non-Javadoc)
	 * @see fi.csc.microarray.messaging.MessagingEndpointIntrfc#replyToMessage(fi.csc.microarray.messaging.message.ChipsterMessage, fi.csc.microarray.messaging.message.ChipsterMessage, java.lang.String)
	 */
    @Override
	public void replyToMessage(ChipsterMessage original, ChipsterMessage reply, String replyChannel) throws JMSException {
    	reply.setMultiplexChannel(replyChannel);
    	Destination replyToDest = original.getReplyTo();
    	sendMessage(replyToDest, reply);
    }


	/**
	 * Not multithread safe.
	 */
    public void sendMessageToClientReplyChannel(Destination replyToDest, ChipsterMessage message) throws JMSException {
		message.setMultiplexChannel(DEFAULT_REPLY_CHANNEL);
		sendMessage(replyToDest, message);
    }

    
    /**
	 * Not multithread safe.
	 */
    private void sendMessage(Destination replyToDest, ChipsterMessage message) throws JMSException {
		Session session = connection.createSession(false, Session.AUTO_ACKNOWLEDGE);
    	try {
			MapMessage mapMessage = session.createMapMessage();
	    	message.marshal(mapMessage);
	    	session.createProducer(replyToDest).send(mapMessage);
    	} finally {
    		session.close();
    	}
    }
    
	/**
	 * Reacts to administrative messages.
	 */
	public void onChipsterMessage(ChipsterMessage msg) {
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

	/* (non-Javadoc)
	 * @see fi.csc.microarray.messaging.MessagingEndpointIntrfc#close()
	 */
    @Override
	public void close() throws JMSException {
    	connection.stop();
    	connection.close();	
    }
    
	/* (non-Javadoc)
	 * @see fi.csc.microarray.messaging.MessagingEndpointIntrfc#getAuthenticationListener()
	 */
	@Override
	public AuthenticationRequestListener getAuthenticationListener() {
		return authenticationListener;
	}
	
	/* (non-Javadoc)
	 * @see fi.csc.microarray.messaging.MessagingEndpointIntrfc#setAuthenticationListener(fi.csc.microarray.messaging.auth.AuthenticationRequestListener)
	 */
	@Override
	public void setAuthenticationListener(AuthenticationRequestListener authenticationListener) {
		this.authenticationListener = authenticationListener;
	}

	/* (non-Javadoc)
	 * @see fi.csc.microarray.messaging.MessagingEndpointIntrfc#getSessionID()
	 */
	@Override
	public String getSessionID() {
		return sessionID;
	}

	/* (non-Javadoc)
	 * @see fi.csc.microarray.messaging.MessagingEndpointIntrfc#setSessionID(java.lang.String)
	 */
	@Override
	public void setSessionID(String sessionID) {
		this.sessionID = sessionID;
	}

}
