package fi.csc.microarray.analyser;

import java.io.IOException;

import javax.jms.JMSException;

import fi.csc.microarray.MicroarrayConfiguration;
import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.util.config.ConfigurationLoader.OldConfigurationFormatException;

/**
 * Main class for launching AnalyserServer.
 * 
 *  @author hupponen
 *
 */
public class AnalyserMain {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws MicroarrayException 
	 * @throws JMSException 
	 * @throws OldConfigurationFormatException 
	 */
	public static void main(String[] args) throws IOException, JMSException, MicroarrayException, OldConfigurationFormatException {
		MicroarrayConfiguration.loadConfiguration();
		new AnalyserServer();	
	}

}
