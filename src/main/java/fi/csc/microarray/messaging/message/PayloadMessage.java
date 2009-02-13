/*
 * Created on Feb 11, 2005
 *
 *
 */
package fi.csc.microarray.messaging.message;

import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.UUID;

import javax.jms.JMSException;
import javax.jms.MapMessage;

import org.apache.log4j.Logger;

import fi.csc.microarray.MicroarrayConfiguration;
import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.filebroker.FileBrokerConfig;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.UrlTransferUtil;
import fi.csc.microarray.util.IOUtils.CopyProgressListener;


/**
 * For sending jobs to back-end components.
 * 
 * @author Taavi Hupponen, Aleksi Kallio
 *
 */
public class PayloadMessage extends ParameterMessage {

	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(PayloadMessage.class);

	private static final String KEY_PAYLOAD_PREFIX = "payload_";
	
	private Map<String, URL> payloads = new HashMap<String, URL>();

	private static boolean useChunked;
	
	static {
		 String chunkedConfig = MicroarrayConfiguration.getValue("messaging", "use_chunked_http");
		 // use chunked if not explicitly disabled
		 useChunked = ! (chunkedConfig != null && chunkedConfig.equals("false")); 
	}
	
	/**
	 * For reflection compatibility (newInstance). DO NOT REMOVE!
	 */
	public PayloadMessage() {
		super();
	}
	
	public PayloadMessage(List<String> parameters) {
		super(parameters);
	}

	@Override
	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);

		// load payload urls
		try {
			for (Enumeration names = from.getMapNames(); names.hasMoreElements(); ) {
				String name = (String)names.nextElement();
				logger.debug("examining " + name);
				if (name.startsWith(KEY_PAYLOAD_PREFIX)) {
					String payloadName = name.substring(KEY_PAYLOAD_PREFIX.length());
					URL url = new URL(from.getString(name));
					payloads.put(payloadName, url);
					logger.debug("Unmarshalled " + name + " -> " + payloadName + ", " + url.toExternalForm());
				}
			}
		} catch (Exception e) {
			handleException(e);
		}
	}

	
	/**
	 * Construct a MapMessage that can be used to create a new JobMessage.
	 */
	@Override
	public void marshal(MapMessage mapMessage) throws JMSException {
		super.marshal(mapMessage);
		
		// add payload urls
		String key;
		String urlString;
		try {
			for (String name : payloadNames()) {
				key = KEY_PAYLOAD_PREFIX + name;
				urlString = payloads.get(name).toExternalForm();
				mapMessage.setString(key, urlString);
				logger.debug("Marshalled " + name + " -> " + key + " " + urlString);
			}
		} catch (Exception e) {
			handleException(e);
		}
		
	}

	
	/**
	 * Returns payload for the job.
	 */
	public InputStream getPayload(String name) throws JMSException {
		if (payloads.get(name) != null) {

			URL url = payloads.get(name);
			InputStream payload = null;
			logger.debug("Opening input stream for " + name + " " + url.toExternalForm());
			try {
				payload = url.openStream();
				
				// make sure that the payload actually is available, as
				// jetty is sometimes a bit slow to write the payload to disk
				// on server side this isn't actually useful as 
				// OnDiskAnalysisJobBase checks that content is available
				int maxWaitTime = 30000;
				int waitTime = 10;
				while (payload.available() <= 0 && waitTime < maxWaitTime) {
					logger.debug("Waiting for payload to become available.");
					waitTime = waitTime*2;
					try {
						Thread.sleep(waitTime);
					} catch (InterruptedException e) {
						logger.error("Interrupted while waiting for payload to become available.");
					}
				}
				
				// rather return no payload than empty payload
				if (payload.available() <= 0) {
					payload = null;
				}
				
			} catch (IOException ioe) {
				logger.error("connection refused when fetching data from URL " + url);
				handleException(ioe);
			}
			return payload;
		
		} else {
			return null;
		}
	}

	
	
	public URL getPayloadURL(String name) {
		return payloads.get(name);
	}
	
	
	/**
	 * Add a new payload to the job message. In practice, upload the payload to server side and
	 * add the URL pointing to it to the list payloads of the message.
	 * 
	 * Payload is uploaded every time, no caching is used.
	 * 
	 * 
	 * @param name name of the payload (the input name of the operation, not the name of the databean)
	 * @param payload
	 * @return the URL pointing to the uploaded payload on the server side
	 * @throws JMSException
	 */
	public URL addPayload(String name, InputStream payload, CopyProgressListener progressListener) throws JMSException {
		URL url = null;

		logger.debug("adding payload " + name);
		// create url and upload the payload to fileserver
		try {
			url = createPayloadUrl(name);
			UrlTransferUtil.uploadStream(url, payload, useChunked, progressListener); // uploadStream flushes and closes
			logger.debug("successfully uploaded " + name + " to " + url);

		} catch (Exception e) {
			handleException(e);
		}
		
		// add the url of uploaded file
		payloads.put(name, url);
		
		return url;
	}
	
	/**
	 * Add the URL of an already uploaded payload to the payloads.
	 * 
	 * Used when a cached copy of the payload already exists on the server side. 
	 * 
	 * @param payloadName name of the payload (the input name of the operation, not the name of the databean)
	 * @param payloadURL an URL pointing to the server side location of the payload
	 */
	public void addPayload(String payloadName, URL payloadURL) {
		payloads.put(payloadName, payloadURL);
	}
	

	/**
	 * Add the contents of a DataBean as a payload. Upload the contents to the server side only
	 * if it is not there already.
	 * 
	 * Control is passed back to the DataBean, takes care of synchronisation issues and also
	 * knows if the contents of the DataBean have been changed.
	 * 
	 * All the messaging stuff is still handled by PayloadMessage, which is passed along to 
	 * the DataBean.
	 * 
	 * 
	 * @param payloadName name of the payload (the input name of the operation, not the name of the databean)
	 * @param dataBean
	 * @throws JMSException
	 * @throws MicroarrayException
	 * @throws IOException
	 */
	public void addPayload(String payloadName, DataBean dataBean, IOUtils.CopyProgressListener progressListener) throws JMSException, MicroarrayException, IOException {
		try {
			dataBean.updateRemoteCache(payloadName, this, progressListener);
		} catch (Exception e) {
			handleException(e);
		}
	}	
	
	
	
	public Set<String> payloadNames() {
		return payloads.keySet();
	}
	
	@Override
	public String toString() {
		return super.toString() + ", payload count: " + payloadNames().size();		
	}

	public Map<String, URL> getPayloads() {
		return payloads;
	}

	public void setPayloads(Map<String, URL> payloads) {
		this.payloads = payloads;
	}

	
	/**
	 * Check if the server side copy of the payload is ok.
	 * 
	 * For now, check that payload exists and that the content length matches
	 * that of the local content.
	 * 
	 * 
	 * @param cachedURL
	 * @param contentLength
	 * @return
	 */
	public boolean checkCachedPayload(URL cachedURL, long contentLength) {
		
		HttpURLConnection connection = null;
		try {

			connection = (HttpURLConnection) cachedURL.openConnection();
			
			// check file existence
			if (connection.getResponseCode() != HttpURLConnection.HTTP_OK) {
				return false;
			} 

		} catch (IOException ioe) {
			return false;
			
		} finally {
			IOUtils.disconnectIfPossible(connection);
		}
		return true;
	}


	/**
	 * Create an unique URL for server side copy of a payload.
	 * 
	 * @param name of the payload
	 * @return
	 * @throws MalformedURLException
	 */
	private static URL createPayloadUrl(String name) throws MalformedURLException {
		String filename = UUID.randomUUID().toString() + "-" + name;
		return FileBrokerConfig.generateUrlToSomeFileBroker(filename);
	}


}
