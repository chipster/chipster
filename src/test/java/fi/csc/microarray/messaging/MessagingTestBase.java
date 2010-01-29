package fi.csc.microarray.messaging;

import org.apache.log4j.Logger;
import org.testng.annotations.AfterSuite;
import org.testng.annotations.BeforeSuite;

import fi.csc.microarray.DemoAuthenticationRequestListener;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.messaging.auth.SimpleAuthenticationRequestListener;

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

	
	
	@BeforeSuite
	protected void setUp() throws Exception {
		
		// use demo listener if no username or password
		System.out.println("setting up authentication listener");
		if (this.username == null || this.password == null) {
			this.authenticationListener = new DemoAuthenticationRequestListener();
		} else {
			this.authenticationListener = new SimpleAuthenticationRequestListener(this.username, this.password);
		}
		
		// use default config if no config url
		System.out.println("initialising client");
		if (configURL == null) {
			DirectoryLayout.initialiseClientLayout().getConfiguration();
		} else {
			DirectoryLayout.initialiseClientLayout(configURL);
		}
		
		System.out.println("setting up messaging");
		logger = Logger.getLogger(MessagingTestBase.class);
		logger.debug("loaded config");
		endpoint =  new MessagingEndpoint(new NodeBase() {
			public String getName() {
				return "test";
			}
		});
		logger.debug("endpoint created");
		endpoint.setAuthenticationListener(authenticationListener);
	}
	
	@AfterSuite
	protected void tearDown() throws Exception {
		System.out.println("closing messaging endpoint");
		endpoint.close();
		endpoint = null;
	}

}
