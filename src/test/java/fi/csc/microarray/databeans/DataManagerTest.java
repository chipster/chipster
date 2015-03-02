package fi.csc.microarray.databeans;

import java.io.File;
import java.io.IOException;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import fi.csc.microarray.ClientContextUtil;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.session.SessionManager;

public class DataManagerTest {

	private DataManager manager;
	private SessionManager sessionManager; 
	
	@Before
	public void init() throws Exception {
		ClientContextUtil.setupClientContext();
		this.manager = Session.getSession().getDataManager();
		this.sessionManager = new SessionManager(manager, null, Session.getSession().getServiceAccessor().getFileBrokerClient(), null);
	}
	
	@Test
	public void testDataTypes() throws IOException {
		File file = File.createTempFile("test", ".png");
		Assert.assertEquals(manager.guessContentType(file).getType(), "image/png");
	}

	
	@Test
	public void testRemoteSessions() throws Exception {
		
		// populate with crap
		File content = File.createTempFile("content", ".tsv");
		DataBean data = manager.createDataBean("test-content", content);
		ClientContextUtil.setupDatabean(data);
		manager.connectChild(data, manager.getRootFolder());
		
		// save
		File session = File.createTempFile("test-remote-session", ".zip");
		sessionManager.saveLightweightSession(session);

		// clear
		manager.deleteAllDataItems();

		// load
		sessionManager.loadLocalSession(session, true);
		
		// check
		Assert.assertEquals(manager.getRootFolder().getChildCount(), 1);
		Assert.assertEquals(manager.getRootFolder().getChildren().iterator().next().getName(), "test-content");
		
		// clean up
		manager.deleteAllDataItems();
	}

}
