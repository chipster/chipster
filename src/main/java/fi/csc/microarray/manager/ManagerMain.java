package fi.csc.microarray.manager;

import java.io.IOException;
import java.sql.SQLException;

import javax.jms.JMSException;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.config.ConfigurationLoader.OldConfigurationFormatException;

/**
 * Main class for launching Manager.
 * 
 * Needed because loading the configuration does not work very
 * well if done in the Manager class (logger problems).
 * 
 *  @author Taavi Hupponen
 *
 */
public class ManagerMain {

	public static void main(String[] args) throws IOException, JMSException, MicroarrayException, OldConfigurationFormatException, ClassNotFoundException, SQLException {
		DirectoryLayout.initialiseServerLayout().getConfiguration();			
		new Manager();	
	}
}
