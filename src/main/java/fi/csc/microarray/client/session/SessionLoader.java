package fi.csc.microarray.client.session;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.zip.ZipFile;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataBean.StorageMethod;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.XmlUtil;

public class SessionLoader {
	
	private File sessionFile;
	private Document sessionDoc;
	private Element sessionElement;
	
	private LinkedHashMap<String, DataFolder> folders = new LinkedHashMap<String, DataFolder>();
	private HashMap<DataFolder, Element> folderElements = new HashMap<DataFolder, Element>();

	private LinkedHashMap<String, DataBean> dataBeans = new LinkedHashMap<String, DataBean>();
	private HashMap<DataBean, Element> dataBeanElements = new HashMap<DataBean, Element>();

	private LinkedHashMap<String, Operation> operations = new LinkedHashMap<String, Operation>();
	private HashMap<Operation, Element> operationElements = new HashMap<Operation, Element>();

	
	private DataManager dataManager;
	
	private static final Logger logger = Logger.getLogger(SessionLoader.class);
	
	
	public SessionLoader(File sessionFile) throws MicroarrayException {
		if (!ClientSession.isValidSessionFile(sessionFile)) {
			throw new MicroarrayException("Not a valid session file.");
		}
		this.sessionFile = sessionFile;
		
		this.dataManager = Session.getSession().getDataManager();
	}
	
	/**
	 * For testing.
	 * 
	 * @param metadataStream
	 * @throws ParserConfigurationException 
	 * @throws IOException 
	 * @throws SAXException 
	 */
	SessionLoader(InputStream metadataStream) throws SAXException, IOException, ParserConfigurationException {
		// create the dom for session.xml
		this.sessionDoc = XmlUtil.parseReader(new BufferedReader(new InputStreamReader(metadataStream)));
		this.sessionElement = sessionDoc.getDocumentElement();
	
		this.dataManager = new DataManager();
	}
	
	public void loadSession() {
		ZipFile zipFile = null;
		try {
			// get the session.xml zip entry
			zipFile = new ZipFile(sessionFile);
			Reader metadataReader = new BufferedReader(new InputStreamReader(zipFile.getInputStream(zipFile.getEntry(ClientSession.SESSION_METADATA_FILENAME))));

			// create the dom for session.xml
			this.sessionDoc = XmlUtil.parseReader(metadataReader);
			this.sessionElement = sessionDoc.getDocumentElement();
			
			
			parseFolders();
			
			for (String folderId : this.folders.keySet()) {
				dataManager.getRootFolder().addChild(folders.get(folderId));
			}
		} 
		// FIXME
		catch (Exception e) {
			logger.error(e);
		}

		// try to close all input streams from the zip file
		finally {
			if (zipFile != null) {
				try {
					zipFile.close();
				} catch (IOException e) {
					logger.warn("could not close zip file");
				}
			}
		}
	}

	void parseFolders() {
		for (Element folderElement : XmlUtil.getChildElements(sessionElement, ClientSession.ELEMENT_FOLDER)) {
			String name = XmlUtil.getChildElement(folderElement, ClientSession.ELEMENT_NAME).getTextContent();
			String id = XmlUtil.getChildElement(folderElement, ClientSession.ELEMENT_ID).getTextContent();
			
			// check for unique id
			if (folders.containsKey(id)) {
				logger.warn("duplicate folder id: " + id + " , ignoring folder: " + name);
				continue;
			}
			
			// create the folder
			DataFolder dataFolder = dataManager.createFolder(name);
			folders.put(id, dataFolder);
			folderElements.put(dataFolder, folderElement);
	
			logger.debug("successfully parsed folder element: " + dataFolder.getName());
		}
	}

	void parseDataBeans() {
		for (Element element : XmlUtil.getChildElements(sessionElement, "data")) {
			String name = XmlUtil.getChildElement(element, ClientSession.ELEMENT_NAME).getTextContent();
			String id = XmlUtil.getChildElement(element, ClientSession.ELEMENT_ID).getTextContent();
			
			// check for unique id
			if (dataBeans.containsKey(id)) {
				logger.warn("duplicate data bean id: " + id + " , ignoring data bean: " + name);
				continue;
			}
			
			// create the data bean
			String storageMethodString = XmlUtil.getChildElement(element, ClientSession.ELEMENT_STORAGE_METHOD).getTextContent();
			StorageMethod storageMethod = DataBean.StorageMethod.valueOf(storageMethodString);
			String urlString = XmlUtil.getChildElement(element, ClientSession.ELEMENT_URL).getTextContent();
			URL url = null;
			try {
				url = new URL(urlString);
			} catch (MalformedURLException e) {
				logger.warn("could not parse url: "  + urlString + " for data bean: " + name);
			}
			
			String cacheURLString = XmlUtil.getChildElement(element, ClientSession.ELEMENT_CACHE_URL).getTextContent();
			URL cacheURL = null;
			try {
				cacheURL = new URL(cacheURLString);
			} catch (MalformedURLException e1) {
				logger.warn("could not parse cache url: "  + cacheURLString + " for data bean: " + name);
			}
			
			DataBean dataBean = null;
			switch (storageMethod) {
			case LOCAL_SESSION:
				try {
					dataBean = dataManager.createDataBeanFromZip(name, url);
				} catch (MicroarrayException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
				break;
			case LOCAL_USER:
				try {
					dataBean = dataManager.createDataBean(name, url);
				} catch (MicroarrayException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				break;
			default:
				logger.warn("unexpected storage method " + storageMethod.name() + " for data bean: " + name);	
				continue;
			}

			dataBean.setCacheUrl(cacheURL);
			dataBean.setNotes(XmlUtil.getChildElement(element, ClientSession.ELEMENT_NOTES).getTextContent());
			//			dataBean.setCreationDate(date);
			
			dataBeans.put(id, dataBean);
			dataBeanElements.put(dataBean, element);
	
			logger.debug("successfully parsed databean element: " + dataBean.getName());
		}
	}

	
	void parseOperations() {
		for (Element element : XmlUtil.getChildElements(sessionElement, ClientSession.ELEMENT_OPERATION)) {
			String id = XmlUtil.getChildElement(element, ClientSession.ELEMENT_ID).getTextContent();
			
			// check for unique id
			if (operations.containsKey(id)) {
				logger.warn("duplicate operation id: " + id);
				continue;
			}

			Element descriptionElement = XmlUtil.getChildElement(element, ClientSession.ELEMENT_DESCRIPTION);
			Element nameElement = XmlUtil.getChildElement(descriptionElement, ClientSession.ELEMENT_NAME);
			String operationId = XmlUtil.getChildElement(nameElement, ClientSession.ELEMENT_ID).getTextContent();
			String displayName = XmlUtil.getChildElement(nameElement, ClientSession.ELEMENT_DISPLAY_NAME).getTextContent();
			
			
			// create the operation
//			Operation operation = new Operation();
//			folders.put(id, dataFolder);
//			folderElements.put(dataFolder, folderElement);
//	
//			logger.debug("successfully parsed folder element: " + dataFolder.getName());
		}
	}

	

}
