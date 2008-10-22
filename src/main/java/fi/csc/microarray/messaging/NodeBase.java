/*
 * Created on Jan 28, 2005
 *
 */
package fi.csc.microarray.messaging;

import org.apache.log4j.Logger;

import java.net.InetAddress;
import java.net.UnknownHostException;

import javax.jms.JMSException;

/**
 * Base class for implementing Node. This class is thread safe.
 * 
 * @author akallio
 */
public abstract class NodeBase implements Node {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(NodeBase.class);

	/**
	 * Uses InetAddress to get the hostname.
	 */
	public String getHost() {
		try {
			return InetAddress.getLocalHost().getHostName();
		} catch (UnknownHostException e) {
			return "unknown";
		}
	}
	
	public void onException(JMSException e) {
		logger.error(e, e);
	}


}
