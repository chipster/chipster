package fi.csc.microarray.messaging;

import java.util.Arrays;

import org.apache.log4j.Logger;
import org.junit.After;
import org.junit.Before;

import fi.csc.microarray.DemoAuthenticationRequestListener;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.messaging.auth.SimpleAuthenticationRequestListener;
import fi.csc.microarray.util.Exceptions;

public abstract class MessagingTestBase {
	/**
	 * Logger for this class
	 */
	private static Logger logger;
	
	protected MessagingEndpoint endpoint;
	protected SimpleAuthenticationRequestListener authenticationListener;
	private String configURL;
	private String username;
	private String password;
	
	public MessagingTestBase() {
	}
	
	public MessagingTestBase(String username, String password) {
		this.username = username;
		this.password = password;
	}

	public MessagingTestBase(String username, String password, String configURL) {
		this.configURL = configURL;
		this.username = username;
		this.password = password;
	}

	
	
	@Before
	public void setUp() throws Exception {
		
		
		// use demo listener if no username or password
		System.out.println("setting up authentication listener");
		if (this.username == null || this.password == null) {
			this.authenticationListener = new DemoAuthenticationRequestListener();
		} else {
			this.authenticationListener = new SimpleAuthenticationRequestListener(this.username, this.password);
		}
		
		// use default config if no config url
		System.out.println("initialising client");
		DirectoryLayout.uninitialise();
		if (configURL == null) {
			DirectoryLayout.initialiseServerLayout(Arrays.asList(new String[]{"client"}));
		} else {
			try {
				DirectoryLayout.initialiseClientLayout(configURL);
			} catch (IllegalStateException e) {
				// ignore
			}
		}
		
		System.out.println("setting up messaging");
		logger = Logger.getLogger(MessagingTestBase.class);
		logger.debug("loaded config");
		endpoint =  new JMSMessagingEndpoint(new NodeBase() {
			public String getName() {
				return "test";
			}
		});
		System.out.println("endpoint created");
		endpoint.setAuthenticationListener(authenticationListener);		
	}
	
	@After
	public void tearDown() {
		if (endpoint != null) {
			System.out.println("closing messaging endpoint");
			try {
				endpoint.close();
			} catch (Exception e) {
				System.out.println("got exception when closing messaging endpoint");
				System.out.println(Exceptions.getStackTrace(e));
			}
			endpoint = null;
		}
	}
}
