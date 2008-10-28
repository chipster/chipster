package fi.csc.microarray.messaging;

import org.apache.log4j.Logger;
import org.testng.annotations.AfterSuite;
import org.testng.annotations.BeforeSuite;

import fi.csc.microarray.DemoAuthenticationRequestListener;
import fi.csc.microarray.MicroarrayConfiguration;

public abstract class MessagingTestBase {
	/**
	 * Logger for this class
	 */
	private static Logger logger;
	
	protected MessagingEndpoint endpoint;
	protected DemoAuthenticationRequestListener authenticationListener = new DemoAuthenticationRequestListener();
	
	@BeforeSuite
	protected void setUp() throws Exception {
		MicroarrayConfiguration.loadConfiguration();
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
