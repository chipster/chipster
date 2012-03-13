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
	private static Logger logger = null;

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
		// lazy init to prevent trouble with unitialised working dirs
		if (logger == null) {
			logger = Logger.getLogger(NodeBase.class);
		}
		// Don't log these as they might also be ok, for example when sending
		// to a deleted temp topic
		// logger.error(e, e);
	}


}
