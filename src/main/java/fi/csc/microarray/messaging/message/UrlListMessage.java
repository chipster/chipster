package fi.csc.microarray.messaging.message;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;

import javax.jms.JMSException;
import javax.jms.MapMessage;


public class UrlListMessage extends ChipsterMessage {
	
	private final static String KEY_URL_LIST = "url-list";

	//this should be safe because other line feeds are encoded to %0A
	private static final String URL_DELIMITER = "\n";
	
	private List<URL> urlList;
	
	public UrlListMessage() {
		super();
	}
	
	/**
	 * Create a new message with the supplied url list. Url items must be in <u>encoded</u> form. 
	 * 
	 * @param url
	 */
	public UrlListMessage(List<URL> url) {
		super();
		this.urlList = url;
	}
	
	public List<URL> getUrlList() {
		return urlList;
	}
	
	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);
		try {
			
			String messageString = from.getString(KEY_URL_LIST);
			String[] urlArray = messageString.split(URL_DELIMITER);
			
			urlList = new LinkedList<URL>();
			
			for (String urlString : urlArray) {			
				urlList.add(new URL(urlString));
			}
			
		} catch (MalformedURLException e) {
			handleException(e);
		}
	}

	public void marshal(MapMessage mapMessage) throws JMSException {
		super.marshal(mapMessage);
		
		String messageString = "";
		
		for (URL url : urlList) {
			messageString += (url.toString() + URL_DELIMITER);
		}
		
		mapMessage.setString(KEY_URL_LIST, messageString);
	}
}
