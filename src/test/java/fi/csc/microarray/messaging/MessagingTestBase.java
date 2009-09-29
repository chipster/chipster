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
	protected SimpleAuthenticationRequestListener authenticationListener = new DemoAuthenticationRequestListener();
	private String configURL;
	
	public MessagingTestBase() {
		this.authenticationListener = new DemoAuthenticationRequestListener();
		configURL = null;
	}
	
	public MessagingTestBase(String username, String password) {
		this.authenticationListener = new SimpleAuthenticationRequestListener(username, password);
		configURL = null;
	}

	public MessagingTestBase(String username, String password, String configURL) {
		this.authenticationListener = new SimpleAuthenticationRequestListener(username, password);
		this.configURL = configURL;
	}

	@BeforeSuite
	protected void setUp() throws Exception {
		if (configURL == null) {
			DirectoryLayout.initialiseClientLayout().getConfiguration();
		} else {
			DirectoryLayout.initialiseClientLayout(configURL);
		}
		
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
		endpoint.close();
		endpoint = null;
	}

}
