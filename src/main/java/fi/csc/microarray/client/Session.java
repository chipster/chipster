/*
 * Created on Feb 3, 2005
 *
 */
package fi.csc.microarray.client;

import java.util.List;

import fi.csc.microarray.client.visualisation.VisualisationMethodRepository;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.module.Module;
import fi.csc.microarray.module.ModuleManager;



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
	private ModuleManager moduleManager;
	private Frames frames;
	private VisualisationMethodRepository visualisations = new VisualisationMethodRepository();
	private String username;
	
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

	public List<Module> getModules() {
		return moduleManager.getModules();
	}

	public Module getPrimaryModule() {
		return moduleManager.getPrimaryModule();
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

	public void setModuleManager(ModuleManager modules) {
		this.moduleManager = modules;
	}

	public Frames getFrames() {
		return frames;
	}

	public void setFrames(Frames frames) {
		this.frames = frames;
	}

	public VisualisationMethodRepository getVisualisations() {
		return visualisations;
	}

	public String getUsername() {
		return username;
	}

	public void setUsername(String username) {
		this.username = username;
	}
	
}
