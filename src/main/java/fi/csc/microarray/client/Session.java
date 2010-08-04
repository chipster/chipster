/*
 * Created on Feb 3, 2005
 *
 */
package fi.csc.microarray.client;

import java.util.HashMap;
import java.util.Map;

import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.module.Modules;



/**
 * Simple holder for 
 * @author Aleksi Kallio
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

	public ClientApplication getApplication() {
		return (ClientApplication)sessionBoundObjects.get("application");
	}
	
	public DataManager getDataManager() {
		return (DataManager)sessionBoundObjects.get("data-manager");
	}

	public ServiceAccessor getServiceAccessor() {
		return (ServiceAccessor)sessionBoundObjects.get("service-accessor");		
	}

	public Modules getModules() {
		return (Modules)sessionBoundObjects.get("modules");
	}
}
