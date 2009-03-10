package fi.csc.microarray.manager;

import java.io.IOException;
import java.sql.SQLException;

import javax.jms.JMSException;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;

// FIXME obsolete, should be removed
public class ManagerMain {

	public static void main(String[] args) throws IOException, JMSException, MicroarrayException, IllegalConfigurationException, ClassNotFoundException, SQLException {
		new Manager();	
	}
}
