package fi.csc.microarray.filebroker;

import java.io.File;
import java.net.MalformedURLException;
import java.util.regex.Pattern;

import fi.csc.microarray.client.RemoteServiceAccessor;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.session.SessionManager;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.module.ModuleManager;

/**
 * Sessions can be stored as server sessions or zip files. This class 
 * converts a session between those two formats.
 * 
 * The client functionality of loading and saving sessions is utilized here,
 * but the client code communicates directly with filebroker using  DirectFileBrokerEndpoint, instead
 * of JMS mediated communication. This eliminates the need for messaging configuration and authentication.
 * 
 * @author klemela
 */
public class ServerSessionImportExportTool {
	
	private static final String LOCAL_FILEBROKER_IP = "127.0.0.1";

	public DataManager dataManager;

	private SessionManager sessionManager;
	
	//A whitelist of allowed characters in file names 
	private static Pattern nameCheck = Pattern.compile("[a-zA-Z0-9_\\-\\.\\ \\(\\)]*");

	public ServerSessionImportExportTool(MessagingEndpoint endpoint) throws Exception {
			
		//set up some client parts
		this.dataManager = new DataManager();
		//setup service accessor manually to customize messaging endpoint
		RemoteServiceAccessor serviceAccessor = new RemoteServiceAccessor();
		FileBrokerClient fileBrokerClient = new JMSFileBrokerClient(
				endpoint.createTopic(Topics.Name.FILEBROKER_TOPIC, AccessMode.WRITE), 
				null, 
				LOCAL_FILEBROKER_IP );
		serviceAccessor.initialise(endpoint, dataManager, fileBrokerClient);
		
		this.sessionManager = new SessionManager(dataManager, null, fileBrokerClient, null);
		//module manager is needed when session loading checks if a file is phenodata
		ModuleManager moduleManager = new ModuleManager();		
		moduleManager.plugAll(dataManager, null);

		Session.getSession().setDataManager(dataManager);	
		Session.getSession().setServiceAccessor(serviceAccessor);
		Session.getSession().setModuleManager(moduleManager);
	}
		
	/**
	 * Load a zip session and save it as a storage session.
	 * 
	 * @param session
	 * @throws Exception
	 */
	public void importSession(File session) throws Exception {

		this.sessionManager.loadLocalSession(session, false);		
		String saveName = filenameToBasename(session.getName());		
		this.sessionManager.saveStorageSession(saveName);		
		this.dataManager.deleteAllDataItems();
	}

	/**
	 * Load a storage session and save it as a zip file.
	 * 
	 * @param uuid
	 * @param zipFile
	 * @throws MalformedURLException
	 * @throws Exception
	 */
	public void exportSession(String uuid, File zipFile) throws MalformedURLException, Exception {

		this.sessionManager.loadStorageSession(uuid);
		this.sessionManager.saveSession(zipFile);		
		this.dataManager.deleteAllDataItems();			
	}
	
	/**
	 * Convert a name of the server session (without directories) to a name of the zip session. 
	 * 
	 * @param dbSessionName
	 * @return
	 * @throws IllegalArgumentException if the name contains illegal characters
	 */
	public String basenameToFileName(String dbSessionName) throws IllegalArgumentException {
		if (nameCheck.matcher(dbSessionName).matches()) {
			return dbSessionName + ".zip";
		} else {
			throw new IllegalArgumentException("Illegal character in server session name '" + dbSessionName + "'");
		}							
	}
	
	/**
	 * Convert a name of the zip session to a server session name (without directories).
	 * 
	 * @param filename
	 * @return
	 * @throws IllegalArgumentException if the name contains illegal characters
	 */
	public String filenameToBasename(String filename) throws IllegalArgumentException {
		if (nameCheck.matcher(filename).matches()) {
			// remove .zip at the end of the line
			// but keep everything else, e.g. .zip.zip becomes .zip 
			return filename.replaceAll(".zip$", "");
		} else {
			throw new IllegalArgumentException("Illegal character in filename '" + filename + "'");
		}
	}
}