/*
 * Created on Feb 3, 2005
 *
 */
package fi.csc.microarray.client;

import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.module.Modules;



/**
 * Simple global access point for client session related objects.
 *  
 * @author Aleksi Kallio
 *
 */
public class Session {
	
	private static Session instance = new Session();	
	
	private ClientApplication clientApplication;
	private DataManager dataManager;
	private ServiceAccessor serviceAccessor;
	private Modules modules;
	private Frames frames;
	
	
	/**
	 * 
	 * @return the singleton object
	 */
	public static Session getSession() {
		return instance;
	}
	
	private Session() {
	}

	public ClientApplication getApplication() {
		return clientApplication;
	}

	public DataManager getDataManager() {
		return dataManager;
	}

	public ServiceAccessor getServiceAccessor() {
		return serviceAccessor;
	}

	public Modules getModules() {
		return modules;
	}

	public void setClientApplication(ClientApplication clientApplication) {
		this.clientApplication = clientApplication;
	}

	public void setDataManager(DataManager dataManager) {
		this.dataManager = dataManager;
	}

	public void setServiceAccessor(ServiceAccessor serviceAccessor) {
		this.serviceAccessor = serviceAccessor;
	}

	public void setModules(Modules modules) {
		this.modules = modules;
	}

	public Frames getFrames() {
		return frames;
	}

	public void setFrames(Frames frames) {
		this.frames = frames;
	}
}
