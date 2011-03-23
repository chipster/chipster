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

import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.stream.StreamSource;

import org.apache.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.session.schema.ChildrenType;
import fi.csc.microarray.client.session.schema.FolderType;
import fi.csc.microarray.client.session.schema.ObjectFactory;
import fi.csc.microarray.client.session.schema.SessionType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataBean.StorageMethod;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.XmlUtil;

public class SessionLoader {
	
	private File sessionFile;
	private SessionType sessionType;
	
//	private LinkedHashMap<String, DataFolder> folders = new LinkedHashMap<String, DataFolder>();
	private HashMap<DataFolder, FolderType> folderTypes = new HashMap<DataFolder, FolderType>();

	private LinkedHashMap<String, DataItem> dataItems = new LinkedHashMap<String, DataItem>();

	
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
	 * @throws JAXBException 
	 * @throws SAXException 
	 */
	SessionLoader(InputStream metadataStream) throws IOException, ParserConfigurationException, JAXBException {

		Unmarshaller unmarshaller = ClientSession.getJAXBContext().createUnmarshaller();
		this.sessionType = unmarshaller.unmarshal(new StreamSource(metadataStream), SessionType.class).getValue();
	
		this.dataManager = new DataManager();
	}
	
	public void loadSession() {
		ZipFile zipFile = null;
		try {
			// get the session.xml zip entry
			zipFile = new ZipFile(sessionFile);
			InputStream metadataStream = zipFile.getInputStream(zipFile.getEntry(ClientSession.SESSION_METADATA_FILENAME));

			// create the dom for session.xml
			Unmarshaller unmarshaller = ClientSession.getJAXBContext().createUnmarshaller();
			this.sessionType = unmarshaller.unmarshal(new StreamSource(metadataStream), SessionType.class).getValue();
		
			parseFolders();
			linkChildren(dataManager.getRootFolder());
			
		} 
		// FIXME
		catch (Exception e) {
			e.printStackTrace();
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
		for (FolderType folderType : sessionType.getFolder()) {
			String name = folderType.getName();
			String id = folderType.getId();
			
			// check for unique id
			if (dataItems.containsKey(id)) {
				logger.warn("duplicate folder id: " + id + " , ignoring folder: " + name);
				continue;
			}
			
			// create the folder
			DataFolder dataFolder;
			if (ClientSession.ROOT_FOLDER_ID.equals(id)) {
				dataFolder = dataManager.getRootFolder();
			} else {
				dataFolder = dataManager.createFolder(name);
			}
			dataItems.put(id, dataFolder);
			folderTypes.put(dataFolder, folderType);
	
			logger.debug("successfully parsed folder element: " + dataFolder.getName());
		}
	}

	void parseDataBeans() {
//		for (Element element : XmlUtil.getChildElements(sessionElement, "data")) {
//			String name = XmlUtil.getChildElement(element, ClientSession.ELEMENT_NAME).getTextContent();
//			String id = XmlUtil.getChildElement(element, ClientSession.ELEMENT_ID).getTextContent();
//			
//			// check for unique id
//			if (dataBeans.containsKey(id)) {
//				logger.warn("duplicate data bean id: " + id + " , ignoring data bean: " + name);
//				continue;
//			}
//			
//			// create the data bean
//			String storageMethodString = XmlUtil.getChildElement(element, ClientSession.ELEMENT_STORAGE_METHOD).getTextContent();
//			StorageMethod storageMethod = DataBean.StorageMethod.valueOf(storageMethodString);
//			String urlString = XmlUtil.getChildElement(element, ClientSession.ELEMENT_URL).getTextContent();
//			URL url = null;
//			try {
//				url = new URL(urlString);
//			} catch (MalformedURLException e) {
//				logger.warn("could not parse url: "  + urlString + " for data bean: " + name);
//			}
//			
//			String cacheURLString = XmlUtil.getChildElement(element, ClientSession.ELEMENT_CACHE_URL).getTextContent();
//			URL cacheURL = null;
//			try {
//				cacheURL = new URL(cacheURLString);
//			} catch (MalformedURLException e1) {
//				logger.warn("could not parse cache url: "  + cacheURLString + " for data bean: " + name);
//			}
//			
//			DataBean dataBean = null;
//			switch (storageMethod) {
//			case LOCAL_SESSION:
//				try {
//					dataBean = dataManager.createDataBeanFromZip(name, url);
//				} catch (MicroarrayException e1) {
//					// TODO Auto-generated catch block
//					e1.printStackTrace();
//				}
//				break;
//			case LOCAL_USER:
//				try {
//					dataBean = dataManager.createDataBean(name, url);
//				} catch (MicroarrayException e) {
//					// TODO Auto-generated catch block
//					e.printStackTrace();
//				}
//				break;
//			default:
//				logger.warn("unexpected storage method " + storageMethod.name() + " for data bean: " + name);	
//				continue;
//			}
//
//			dataBean.setCacheUrl(cacheURL);
//			dataBean.setNotes(XmlUtil.getChildElement(element, ClientSession.ELEMENT_NOTES).getTextContent());
//			//			dataBean.setCreationDate(date);
//			
//			dataBeans.put(id, dataBean);
//			dataBeanElements.put(dataBean, element);
//	
//			logger.debug("successfully parsed databean element: " + dataBean.getName());
//		}
	}

	
	void parseOperations() {
//		for (Element element : XmlUtil.getChildElements(sessionElement, ClientSession.ELEMENT_OPERATION)) {
//			String id = XmlUtil.getChildElement(element, ClientSession.ELEMENT_ID).getTextContent();
//			
//			// check for unique id
//			if (operations.containsKey(id)) {
//				logger.warn("duplicate operation id: " + id);
//				continue;
//			}
//
//			Element descriptionElement = XmlUtil.getChildElement(element, ClientSession.ELEMENT_DESCRIPTION);
//			Element nameElement = XmlUtil.getChildElement(descriptionElement, ClientSession.ELEMENT_NAME);
//			String operationId = XmlUtil.getChildElement(nameElement, ClientSession.ELEMENT_ID).getTextContent();
//			String displayName = XmlUtil.getChildElement(nameElement, ClientSession.ELEMENT_DISPLAY_NAME).getTextContent();
//			
//			
//			// create the operation
////			Operation operation = new Operation();
////			folders.put(id, dataFolder);
////			folderElements.put(dataFolder, folderElement);
////	
////			logger.debug("successfully parsed folder element: " + dataFolder.getName());
//		}
	}

	
	private void linkChildren(DataFolder parent) {
		ChildrenType childrenType = folderTypes.get(parent).getChildren();
		
		// no children at all? 
		if (childrenType == null) {
			return;
		}
		
		for (String childId : childrenType.getChild()) {
			
			// check that the referenced data item exists
			DataItem child = dataItems.get(childId);
			if (child == null) {
				logger.warn("child with id: " + childId + " does not exist");
				continue;
			}

			// add as a child
			parent.addChild(child);
			
			// recursively go inside folders
			if (child instanceof DataFolder) {
				linkChildren((DataFolder) child);
			}
		}
	}
	

}
