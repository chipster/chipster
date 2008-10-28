package fi.csc.microarray.messaging;

import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;

/**
 *
 * @author  Aleksi Kallio
 */
public class JmsTestBase extends MessagingTestBase {
    
    
	@BeforeTest
	protected void setUp() throws Exception {
		super.setUp();
	}
	
	@AfterTest
	protected void tearDown() throws Exception {
		super.tearDown();
	}
	
}
