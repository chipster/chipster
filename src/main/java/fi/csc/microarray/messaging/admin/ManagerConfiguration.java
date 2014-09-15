package fi.csc.microarray.messaging.admin;

import java.io.IOException;
import java.util.Arrays;
import java.util.concurrent.locks.ReentrantLock;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;

public class ManagerConfiguration {
	
	static ReentrantLock mutex = new ReentrantLock();

	public static Configuration getConfiguration() throws IOException, IllegalConfigurationException {
		
		// prevent parallel initialization attempts
		mutex.lock();
		try {
			if (!DirectoryLayout.isInitialised()) {
				DirectoryLayout.initialiseServerLayout(Arrays.asList(new String[] {"manager"}));
			}
		} finally {
			mutex.unlock();
		}

		return DirectoryLayout.getInstance().getConfiguration();
	}
	
	public static void init() throws IOException, IllegalConfigurationException {		
		getConfiguration();
	}
}
