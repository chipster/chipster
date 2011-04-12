package fi.csc.microarray.client.session;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
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
import org.xml.sax.SAXException;

import fi.csc.microarray.client.NameID;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.session.schema.DataType;
import fi.csc.microarray.client.session.schema.FolderType;
import fi.csc.microarray.client.session.schema.InputType;
import fi.csc.microarray.client.session.schema.LinkType;
import fi.csc.microarray.client.session.schema.NameType;
import fi.csc.microarray.client.session.schema.OperationType;
import fi.csc.microarray.client.session.schema.ParameterType;
import fi.csc.microarray.client.session.schema.SessionType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataBean.StorageMethod;
import fi.csc.microarray.exception.MicroarrayException;

public class SessionLoader {
	
	private File sessionFile;
	private SessionType sessionType;
	
	private LinkedHashMap<String, DataFolder> folders = new LinkedHashMap<String, DataFolder>();
	private HashMap<DataFolder, FolderType> folderTypes = new HashMap<DataFolder, FolderType>();

	private LinkedHashMap<String, DataBean> dataBeans = new LinkedHashMap<String, DataBean>();
	private HashMap<DataBean, DataType> dataTypes = new HashMap<DataBean, DataType>();

	private LinkedHashMap<String, OperationRecord> operationRecords = new LinkedHashMap<String, OperationRecord>();
	private HashMap<OperationRecord, OperationType> operationTypes = new HashMap<OperationRecord, OperationType>();

	
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
			InputStream metadataStream = zipFile.getInputStream(zipFile.getEntry(ClientSession.SESSION_DATA_FILENAME));

			// create the dom for session.xml
			Unmarshaller unmarshaller = ClientSession.getJAXBContext().createUnmarshaller();
			this.sessionType = unmarshaller.unmarshal(new StreamSource(metadataStream), SessionType.class).getValue();
		
			parseFolders();
			parseDataBeans();
			parseOperations();
			
			linkOperationsToOutputs();
			linkChildren(dataManager.getRootFolder());

			linkDataBeans();

			linkInputsToOperations();

			
			
			// check
			logger.info("after load checking");
			for (DataBean dataBean : dataBeans.values()) {
				if (dataBean.getOperationRecord() == null) {
					logger.info(dataBean.getName() + " has null operation record");
				} 
			}
			
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
			if (getDataItem(id) != null) {
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
			folders.put(id, dataFolder);
			folderTypes.put(dataFolder, folderType);
	
			logger.debug("successfully parsed folder element: " + dataFolder.getName());
		}
	}

	void parseDataBeans() {
		for (DataType dataType : this.sessionType.getData()) {
			String name = dataType.getName();
			String id = dataType.getId();
			
			// check for unique id
			if (getDataItem(id) != null) {
				logger.warn("duplicate data bean id: " + id + " , ignoring data bean: " + name);
				continue;
			}
			
			// create the data bean
			String storageMethodString = dataType.getStorageType();
			StorageMethod storageMethod = DataBean.StorageMethod.valueOf(storageMethodString);
			String urlString = dataType.getUrl();
			URL url = null;
			try {
				url = new URL(urlString);
			} catch (MalformedURLException e) {
				logger.warn("could not parse url: "  + urlString + " for data bean: " + name);
				continue;
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

			// cache url
			String cacheURLString = dataType.getCacheUrl();
			if (cacheURLString != null) {
				URL cacheURL = null;
				try {
					cacheURL = new URL(cacheURLString);
				} catch (MalformedURLException e1) {
					logger.warn("could not parse cache url: "  + cacheURLString + " for data bean: " + name);
				}
				dataBean.setCacheUrl(cacheURL);
			}

			// notes
			dataBean.setNotes(dataType.getNotes());
			//			dataBean.setCreationDate(date);
			
			dataBean.setContentType(Session.getSession().getDataManager().guessContentType(dataBean.getName()));
			
			dataBeans.put(id, dataBean);
			dataTypes.put(dataBean, dataType);
	
			logger.debug("successfully parsed databean element: " + dataBean.getName());
		}
	}

	
	private void parseOperations() {
		for (OperationType operationType : sessionType.getOperation()) {
			String operationSessionId = operationType.getId();

			// check for unique id
			if (operationRecords.containsKey(operationSessionId)) {
				logger.warn("duplicate operation id: " + operationSessionId);
				continue;
			}

			OperationRecord operationRecord = new OperationRecord();

			// name id
			operationRecord.setNameID(createNameID(operationType.getName()));
			
			// category
			operationRecord.setCategoryName(operationType.getCategory());
			String colorString = operationType.getCategoryColor();
			if (colorString != null) {
				operationRecord.setCategoryColor(Color.decode(colorString));
			}
			
			// parameters
			for (ParameterType parameterType : operationType.getParameter()) {
				operationRecord.addParameter(parameterType.getName().getId(), parameterType.getName().getDisplayName(), parameterType.getName().getDescription(), parameterType.getValue());
			}
			
			// TODO source code
			
			operationRecords.put(operationSessionId, operationRecord);
			operationTypes.put(operationRecord, operationType);
		}
	}

	
	private void linkChildren(DataFolder parent) {
		for (String childId : folderTypes.get(parent).getChild()) {
			
			// check that the referenced data item exists
			DataItem child = getDataItem(childId);
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
	
	/**
	 * Link OperationRecords and DataBeans by adding real input DataBean references
	 * to OperationRecords.
	 * 
	 */
	private void linkInputsToOperations() {
		
		for (OperationRecord operationRecord : operationRecords.values()) {
			
			// FIXME add checks
			// get data bean ids from session data
			for (InputType inputType : operationTypes.get(operationRecord).getInput()) {
				
				// find the data bean
				DataBean inputBean = dataBeans.get(inputType.getData());
				
				// skip phenodata, it is bound automatically
				if (inputBean.queryFeatures("/phenodata/").exists()) {
					continue; 
				}

				// add the reference to the operation record
				operationRecord.addInput(createNameID(inputType.getName()), inputBean);
			}
		}
	}

	
	private NameID createNameID(NameType name) {
		return new NameID(name.getId(), name.getDisplayName(), name.getDescription());
	}

	/**
	 * Add links form DataBeans to the OperationRecord which created the DataBean.
	 * 
	 */
	private void linkOperationsToOutputs() {
		
		for (DataBean dataBean : dataBeans.values()) {
			String operationId = dataTypes.get(dataBean).getResultOf();
			OperationRecord operationRecord = operationRecords.get(operationId);
			dataBean.setOperationRecord(operationRecord);
		}
	}

	
	private void linkDataBeans() {
		for (DataBean dataBean : dataBeans.values()) {
			for (LinkType linkType : dataTypes.get(dataBean).getLink()) {
				dataBean.addLink(Link.valueOf(linkType.getType()), dataBeans.get(linkType.getTarget()));
			}
			
		}
	}
	
	/**
	 * 
	 * @param id
	 * @return null if no DataItem for the id is found
	 */
	private DataItem getDataItem(String id) {
		DataItem dataItem = folders.get(id);
		if (dataItem != null) {
			return dataItem;
		} else { 
			return dataBeans.get(id);
		}
	}
}
