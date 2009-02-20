package fi.csc.microarray.manager;

import java.io.IOException;
import java.sql.SQLException;

import javax.jms.JMSException;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.config.MicroarrayConfiguration;
import fi.csc.microarray.config.ConfigurationLoader.OldConfigurationFormatException;

/**
 * Main class for launching Manager.
 * 
 * Needed because loading the configuration does not work very
 * well if done in the Manager class (logger problems).
 * 
 *  @author hupponen
 *
 */
public class ManagerMain {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws MicroarrayException 
	 * @throws JMSException 
	 * @throws OldConfigurationFormatException 
	 * @throws SQLException 
	 * @throws ClassNotFoundException 
	 */
	public static void main(String[] args) throws IOException, JMSException, MicroarrayException, OldConfigurationFormatException, ClassNotFoundException, SQLException {
		MicroarrayConfiguration.loadConfiguration();
		new Manager();	
	}
}
