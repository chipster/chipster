package fi.csc.microarray.filebroker;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import javax.jms.JMSException;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import fi.csc.microarray.filebroker.FileBrokerClient.FileBrokerArea;
import fi.csc.microarray.messaging.MessagingTestBase;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.security.CryptoKey;
import fi.csc.microarray.util.IOUtils;

public class FileBrokerClientTest extends MessagingTestBase {

	private FileBrokerClient fbc;
	
	@Before
	public void setUp() throws Exception {
		super.setUp();
		fbc = new JMSFileBrokerClient(super.endpoint.createTopic(Topics.Name.FILEBROKER_TOPIC, AccessMode.WRITE));
	}

	@Before
	public void tearDown() {
		super.tearDown();
	}
	
	@Test
	public void testFileBroker() throws FileNotFoundException, FileBrokerException, JMSException, IOException {
		File file = new File("src/test/resources/affy_example.cel");
		
		System.out.println("Adding file");
		String dataId = CryptoKey.generateRandom();
		fbc.addFile(dataId, FileBrokerArea.CACHE, new FileInputStream(file), file.length(), null);

		System.out.println("Checking file");
		// Assert.assertTrue(fbc.checkFile(url, file.length()));
		
		System.out.println("Getting file");
		long outputContentLength = 0;
		BufferedInputStream input = new BufferedInputStream(fbc.getInputStream(dataId));
		while (input.read() != -1) {
			outputContentLength++;
		}
		IOUtils.closeIfPossible(input);
		
		Assert.assertEquals(file.length(), outputContentLength);
		
	}
}
