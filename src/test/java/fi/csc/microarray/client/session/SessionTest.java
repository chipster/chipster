package fi.csc.microarray.client.session;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.net.MalformedURLException;
import java.nio.file.Files;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import org.junit.Assert;
import org.junit.Test;

import fi.csc.microarray.client.RemoteServiceAccessor;
import fi.csc.microarray.client.ServiceAccessor;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.ToolModule;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.filebroker.DbSession;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.filebroker.FileServer;
import fi.csc.microarray.filebroker.JMSFileBrokerClient;
import fi.csc.microarray.filebroker.MockJettyFileServer;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.DirectMessagingEndpoint;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.auth.SimpleAuthenticationRequestListener;
import fi.csc.microarray.module.ModuleManager;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.util.IOUtils;

public class SessionTest {

	private String username = "username";
	private String password = "password";
	private SimpleAuthenticationRequestListener authenticationListener;

	
	@Test
	public void test() throws Exception {

		// initialise and configure
		File workDir = Files.createTempDirectory("remote-session-unit-test-temp").toFile();
		new File(workDir, "conf").mkdir();
		File logDir = new File(workDir, "logs");
		logDir.mkdir();
		new File(workDir, "security").mkdir();
		DirectoryLayout.setBaseDirOverride(workDir);

		File configFile = new File(workDir, "conf" + File.separator + "chipster-config.xml");
		IOUtils.copy(new ByteArrayInputStream(getMockConfig().getBytes()), configFile);

		// boot up file server so that it is connected to mock messaging fabric
		DirectMessagingEndpoint endpoint = new DirectMessagingEndpoint();
		new FileServer(null, endpoint, new MockJettyFileServer());

		// test file broker using JMSFileBrokerClient
		MessagingTopic urlTopic = endpoint.createTopic(Topics.Name.FILEBROKER_TOPIC, AccessMode.WRITE);
		FileBrokerClient fbc = new JMSFileBrokerClient(urlTopic);

		// Test just basic messaging. This used to request for new url, but it's not possible anymore in FileBrokerClient 
		boolean diskSpace = fbc.requestDiskSpace(0);
		
		// check results
		Assert.assertTrue(diskSpace);
		
		// check all log files after execution
		for (File logFile : logDir.listFiles()) {
			ByteArrayOutputStream contents = new ByteArrayOutputStream();
			IOUtils.copy(new FileInputStream(logFile), contents);
			Assert.assertFalse(logFile + " should not contain exception", contents.toString().contains("Exception"));
		}		
	}


	private String getMockConfig() {
		return "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n" + 
				"<configuration content-version=\"3\">\n" + 
				"\n" + 
				"	<configuration-module moduleId=\"messaging\">\n" + 
				"\n" + 
				"		<!-- host of message broker (JMS server ActiveMQ) to connect to -->\n" + 
				"		<entry entryKey=\"broker-host\">\n" + 
				"			<value>not defined</value>\n" + 
				"		</entry>\n" + 
				"\n" + 
				"		<!-- protocol used to connect to message broker -->\n" + 
				"		<entry entryKey=\"broker-protocol\">\n" + 
				"			<value>not defined</value>\n" + 
				"		</entry>\n" + 
				"\n" + 
				"		<!-- port used to connect to message broker -->\n" + 
				"		<entry entryKey=\"broker-port\">\n" + 
				"			<value>-1</value>\n" + 
				"		</entry>\n" + 
				"		\n" + 
				"	</configuration-module>\n" + 
				"\n" + 
				"	<configuration-module moduleId=\"security\">\n" + 
				"\n" + 
				"		<!-- username for authenticating connection to broker -->\n" + 
				"		<entry entryKey=\"username\">\n" + 
				"			<value>not defined</value>\n" + 
				"		</entry>\n" + 
				"\n" + 
				"		<!-- password for authenticating connection to broker -->\n" + 
				"		<entry entryKey=\"password\">\n" + 
				"			<value>not defined</value>\n" + 
				"		</entry>\n" + 
				"	\n" + 
				"	</configuration-module>\n" + 
				"\n" + 
				"	<configuration-module moduleId=\"filebroker\">\n" + 
				"	\n" + 
				"	    <!-- url of this file broker instance -->\n" + 
				"		<entry entryKey=\"url\">\n" + 
				"			<value>http://mockhost</value>\n" + 
				"		</entry>\n" + 
				"		\n" + 
				"		<!-- server port to use in this file broker instance -->\n" + 
				"		<entry entryKey=\"port\">\n" + 
				"			<value>-1</value>			\n" + 
				"        </entry>		\n" + 
				"        \n" + 
				"        \n" + 
				"        <entry entryKey=\"metadata-port\">\n" + 
				"			<value>-1</value>\n" + 
				"       		</entry>\n" + 
				"        \n" + 
				"        \n" + 
				"        <entry entryKey=\"clean-up-trigger-limit-percentage\" type=\"int\" description=\"when disk usage reaches this percentage clean up\">\n" + 
				"			<value>22</value>\n" + 
				"		</entry>\n" + 
				"\n" + 
				"		<entry entryKey=\"clean-up-target-percentage\" type=\"int\" description=\"when cleaning up drop disk usage to this percentage\">\n" + 
				"			<value>20</value>\n" + 
				"		</entry>\n" + 
				"\n" + 
				"\n" + 
				"		<entry entryKey=\"clean-up-minimum-file-age\" type=\"int\" description=\"only clean up files older than this, seconds\">\n" + 
				"			<value>259200</value>\n" + 
				"		</entry>\n" + 
				"\n" + 
				"		<entry entryKey=\"minimum-space-for-accept-upload\" type=\"int\" description=\"when client requests for free space, say no if less than this many megabytes would be available after upload, megabytes\">\n" + 
				"			<value>100</value>\n" + 
				"		</entry>\n" + 
				"        \n" + 
				"	</configuration-module>\n" + 
				"\n" + 
				"</configuration>\n";
	}
	
	
	//@Test
	public void testStorageSessions() throws Exception {
		
		// set up modules
		ModuleManager moduleManager = new ModuleManager("fi.csc.microarray.module.chipster.MicroarrayModule");
		Session.getSession().setModuleManager(moduleManager);
		
		// set up system
		this.authenticationListener = new SimpleAuthenticationRequestListener(username, password);
		DataManager manager = new DataManager();
		moduleManager.plugAll(manager, null);
		LinkedList<ToolModule> toolModules = new LinkedList<ToolModule>();
		ServiceAccessor serviceAccessor = new RemoteServiceAccessor();
		serviceAccessor.initialise(manager, authenticationListener);
		serviceAccessor.fetchDescriptions(new MicroarrayModule());
		Session.getSession().setServiceAccessor(serviceAccessor);
		toolModules.addAll(serviceAccessor.getModules());
		FileBrokerClient fileBrokerClient = serviceAccessor.getFileBrokerClient();
		
		SessionManager sessionManager = new SessionManager(manager, null, fileBrokerClient, null);
		
		// create data
		DataBean data = manager.createDataBean("test");
		data.setParent(manager.getRootFolder());
		int dbCountOrig = manager.databeans().size();
		
		// save storage session
		String sessionName = "unit test session " + new Random().nextInt(10000);
		sessionManager.saveStorageSession(sessionName);
				
		// delete all data
		manager.deleteAllDataItems();
		
		// load session
		List<DbSession> sessions = fileBrokerClient.listRemoteSessions();
		String sessionId = findSession(sessionName, sessions);
		sessionManager.loadStorageSession(sessionId);

		// assert loaded data is ok
		Assert.assertEquals(manager.databeans().size(), dbCountOrig);
		
		// remove remote session
		fileBrokerClient.removeRemoteSession(sessionId);

		// assert removal
		List<DbSession> sessions2 = fileBrokerClient.listRemoteSessions();
		Assert.assertNull(findSession(sessionName, sessions2));

	}

	private String findSession(String sessionName,	List<DbSession> sessions) throws MalformedURLException {
		String sessionId = null;
		for (DbSession session : sessions) {
			if (sessionName.equals(session.getName())) {
				sessionId = session.getDataId();
				break;
			}
		}
		return sessionId;
	}
}
