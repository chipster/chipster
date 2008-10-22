/*
 * Created on Feb 3, 2005
 *
 */
package fi.csc.microarray.client;

import java.util.HashMap;
import java.util.Map;

import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.module.Modules;



/**
 * @author akallio
 *
 */
public class Session {
	private static Session instance = new Session();	
	private Map<String, Object> sessionBoundObjects = new HashMap<String, Object>(); 

	/**
	 * 
	 * @return the singleton object
	 */
	public static Session getSession() {
		return instance;
	}
	
	private Session() {
	}
	
	public void putObject(String name, Object endpoint) {
		sessionBoundObjects.put(name, endpoint);
	}

	public MessagingEndpoint getMessagingEndpoint(String name) {
		return (MessagingEndpoint)sessionBoundObjects.get(name);
	}
		
	public TaskExecutor getJobExecutor(String name) {
		return (TaskExecutor)sessionBoundObjects.get(name);
	}

	public ClientApplication getApplication() {
		return (ClientApplication)sessionBoundObjects.get("application");
	}
	
	public DataManager getDataManager() {
		return (DataManager)sessionBoundObjects.get("data-manager");
	}
	
	public Modules getModules() {
		return (Modules)sessionBoundObjects.get("modules");
	}
}
