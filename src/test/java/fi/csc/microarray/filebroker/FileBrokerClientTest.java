package fi.csc.microarray.filebroker;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;

import javax.jms.JMSException;

import org.testng.Assert;
import org.testng.annotations.AfterSuite;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import fi.csc.microarray.messaging.MessagingTestBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.util.IOUtils;

public class FileBrokerClientTest extends MessagingTestBase {

	private FileBrokerClient fbc;
	
	@BeforeSuite(alwaysRun = true)
	protected void setUp() throws Exception {
		super.setUp();
		fbc = new FileBrokerClient(super.endpoint.createTopic(Topics.Name.URL_TOPIC, AccessMode.WRITE));
	}

	@AfterSuite(alwaysRun = true)
	protected void tearDown() throws Exception {
		super.tearDown();
	}
	
	@Test
	public void testFileBroker() throws FileNotFoundException, FileBrokerException, JMSException, IOException {
		File file = new File("src/test/resources/affy_example.cel");
		
		System.out.println("Adding file");
		URL url = fbc.addFile(new FileInputStream(file), null);

		System.out.println("Checking file");
		Assert.assertTrue(fbc.checkFile(url, file.length()));
		
		System.out.println("Getting file");
		long outputContentLength = 0;
		BufferedInputStream input = new BufferedInputStream(fbc.getFile(url));
		while (input.read() != -1) {
			outputContentLength++;
		}
		IOUtils.closeIfPossible(input);
		
		Assert.assertEquals(file.length(), outputContentLength);
		
	}
}
