package fi.csc.microarray.databeans;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.operation.ToolCategory;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.module.ModuleManager;

public class DataManagerTest {

	private DataManager manager; 
	
	@BeforeClass(alwaysRun = true)
	public void init() throws IOException, IllegalConfigurationException, InstantiationException, IllegalAccessException, ClassNotFoundException {
		DirectoryLayout.initialiseSimpleLayout().getConfiguration();			
		this.manager = new DataManager();
		ModuleManager moduleManager = new ModuleManager("fi.csc.microarray.module.chipster.MicroarrayModule");
		moduleManager.plugAll(this.manager, null);
		Session.getSession().setModuleManager(moduleManager);
	}
	
	@Test(groups = {"unit"} )
	public void testDataTypes() throws IOException {
		File file = File.createTempFile("test", ".png");
		Assert.assertEquals(manager.guessContentType(file).getType(), "image/png");
	}

	
	@Test(groups = {"unit"} )
	public void testRemoteSessions() throws Exception {
		
		// populate with crap
		File content = File.createTempFile("content", ".tsv");
		DataBean data = manager.createDataBean("test-content", content);
		OperationRecord operationRecord = new OperationRecord(new Operation(new OperationDefinition("", "", new ToolCategory(""), "", false), new DataBean[] { data } ));
		data.setOperationRecord(operationRecord);
		operationRecord.setModule("");
		manager.getRootFolder().addChild(data);
		
		// save
		File session = File.createTempFile("test-remote-session", ".zip");
		manager.saveLightweightSession(session);

		// clear
		manager.deleteAllDataItems();

		// load
		manager.loadSession(session, true);
		
		// clean up
		manager.deleteAllDataItems();
	}

}
