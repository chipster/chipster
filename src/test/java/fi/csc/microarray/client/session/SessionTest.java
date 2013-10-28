package fi.csc.microarray.client.session;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.nio.file.Files;
import java.util.LinkedList;
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
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.filebroker.FileServer;
import fi.csc.microarray.messaging.MockMessagingEndpoint;
import fi.csc.microarray.messaging.auth.SimpleAuthenticationRequestListener;
import fi.csc.microarray.module.ModuleManager;
import fi.csc.microarray.module.chipster.MicroarrayModule;

public class SessionTest {

	private String username = "username";
	private String password = "password";
	private SimpleAuthenticationRequestListener authenticationListener;

	
	@Test
	public void test() throws IOException {

		// initialise and configure
		File workDir = Files.createTempDirectory("remote-session-unit-test-temp").toFile();
		new File(workDir, "conf").mkdir();
		new File(workDir, "logs").mkdir();
		new File(workDir, "security").mkdir();
		File configFile = new File(workDir, "conf" + File.separator + "chipster-config.xml");
		DirectoryLayout.setBaseDirOverride(workDir);

		// boot up file server so that it is connected to mock messaging fabric
		MockMessagingEndpoint endpoint = new MockMessagingEndpoint();
		new FileServer(null, endpoint);

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
		
		// create data
		DataBean data = manager.createDataBean("test");
		data.setParent(manager.getRootFolder());
		int dbCountOrig = manager.databeans().size();
		
		// save storage session
		String sessionName = "unit test session " + new Random().nextInt(10000);
		manager.saveStorageSession(sessionName);
				
		// delete all data
		manager.deleteAllDataItems();
		
		// load session
		String[][] sessions = fileBrokerClient.listRemoteSessions();
		URL sessionURL = findSession(sessionName, sessions);
		manager.loadStorageSession(sessionURL);

		// assert loaded data is ok
		Assert.assertEquals(manager.databeans().size(), dbCountOrig);
		
		// remove remote session
		fileBrokerClient.removeRemoteSession(sessionURL);

		// assert removal
		String[][] sessions2 = fileBrokerClient.listRemoteSessions();
		Assert.assertNull(findSession(sessionName, sessions2));

	}

	private URL findSession(String sessionName,	String[][] sessions) throws MalformedURLException {
		URL sessionURL = null;
		for (int i = 0; i < sessions[0].length; i++) {
			if (sessionName.equals(sessions[0][i])) {
				sessionURL = new URL(sessions[1][i]);
				break;
			}
		}
		return sessionURL;
	}

	
}
