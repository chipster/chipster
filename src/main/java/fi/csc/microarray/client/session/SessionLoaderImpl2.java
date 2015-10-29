package fi.csc.microarray.client.session;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipInputStream;

import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.datatype.XMLGregorianCalendar;
import javax.xml.transform.stream.StreamSource;

import org.apache.log4j.Logger;
import org.eclipse.jetty.io.WriterOutputStream;

import de.schlichtherle.truezip.zip.ZipFile;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.NameID;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.operation.OperationRecord.ParameterRecord;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.client.session.schema2.DataType;
import fi.csc.microarray.client.session.schema2.FolderType;
import fi.csc.microarray.client.session.schema2.InputType;
import fi.csc.microarray.client.session.schema2.LinkType;
import fi.csc.microarray.client.session.schema2.LocationType;
import fi.csc.microarray.client.session.schema2.NameType;
import fi.csc.microarray.client.session.schema2.OperationType;
import fi.csc.microarray.client.session.schema2.ParameterType;
import fi.csc.microarray.client.session.schema2.SessionType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataManager.StorageMethod;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.util.IOUtils;

public class SessionLoaderImpl2 {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(SessionLoaderImpl2.class);


	private DataManager dataManager;

	private File sessionFile;
	private String sessionId;
	private SessionType sessionType;
	private boolean isDatalessSession;

	private LinkedHashMap<String, DataFolder> folders = new LinkedHashMap<String, DataFolder>();
	private HashMap<DataFolder, FolderType> folderTypes = new HashMap<DataFolder, FolderType>();

	private LinkedHashMap<String, DataBean> dataBeans = new LinkedHashMap<String, DataBean>();
	private HashMap<DataBean, DataType> dataTypes = new HashMap<DataBean, DataType>();

	private LinkedHashMap<String, OperationRecord> operationRecords = new LinkedHashMap<String, OperationRecord>();
	private HashMap<OperationRecord, OperationType> operationTypes = new HashMap<OperationRecord, OperationType>();

	private ZipFile zipFile;
	private ZipInputStream zipStream;


	private Integer xOffset;


	private String sessionNotes;


	public SessionLoaderImpl2(File sessionFile, DataManager dataManager, boolean isDatalessSession) {
		this.sessionFile = sessionFile;
		this.sessionId = null;
		this.dataManager = dataManager;
		this.isDatalessSession = isDatalessSession;
	}

	public SessionLoaderImpl2(String sessionId, DataManager dataManager, boolean isDatalessSession) {
		this.sessionFile = null;
		this.sessionId = sessionId;
		this.dataManager = dataManager;
		this.isDatalessSession = isDatalessSession;
	}

	/**
	 * 
	 * The open file and/or stream is stored in field, remember to close them.
	 * 
	 * <pre>
	 * try {
	 * 
	 * } finally {
	 * 	IOUtils.closeIfPossible(zipFile);
	 * 	IOUtils.closeIfPossible(zipStream);
	 * }
	 * </pre>
	 * 
	 * @param zipEntry
	 * @return
	 * @throws Exception 
	 */
	private InputStream getStreamOfZipEntry(String zipEntry) throws Exception {

		InputStream stream = null;
		
		if (sessionFile != null) {
			
			if (!sessionFile.exists()) {
				throw new IOException("session file does not exist: " + sessionFile);
			}
			
			// get the zip entry using TrueZip
			zipFile = new ZipFile(sessionFile);
			stream = zipFile.getInputStream(zipEntry);
			
		} else if (sessionId != null) {
			FileBrokerClient fileBrokerClient = Session.getSession().getServiceAccessor().getFileBrokerClient();
			InputStream inputStream = fileBrokerClient.getInputStream(sessionId);
			// get the session.xml zip entry using JDK, we don't need large ZIP support here because URL based sessions have no data
			zipStream = new ZipInputStream(inputStream);
			
			ZipEntry entry;
	        while ((entry = zipStream.getNextEntry()) != null) {
	        	if (zipEntry.equals(entry.getName())) {
	        		stream = zipStream; // zip stream is wound to right entry now, use it
	        		break; 
	        	}
	        }
		}		
		return stream;
	}

	private void parseMetadata() throws Exception {

		try {

			InputStream metadataStream = getStreamOfZipEntry(UserSession.SESSION_DATA_FILENAME);
			
			// validate
			//ClientSession.getSchema().newValidator().validate(new StreamSource(metadataStream));
			if (metadataStream == null) {
				throw new ZipException("session file corrupted, entry " + UserSession.SESSION_DATA_FILENAME + " was missing");
			}
			
			// parse the metadata xml to java objects using jaxb
			Unmarshaller unmarshaller = UserSession.getJAXBContext().createUnmarshaller();
			unmarshaller.setSchema(UserSession.getSchema());
			NonStoppingValidationEventHandler validationEventHandler = new NonStoppingValidationEventHandler();
			unmarshaller.setEventHandler(validationEventHandler);
			this.sessionType = unmarshaller.unmarshal(new StreamSource(metadataStream), SessionType.class).getValue();
			
			if (validationEventHandler.hasEvents()) {
				throw new JAXBException("Invalid session file:\n" + validationEventHandler.getValidationEventsAsString());
			}
		}
		finally {
			IOUtils.closeIfPossible(zipFile);
			IOUtils.closeIfPossible(zipStream);
		}
	}

	private void createFolders() {
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
			if (UserSession.ROOT_FOLDER_ID.equals(id)) {
				dataFolder = dataManager.getRootFolder();
			} else {
				dataFolder = dataManager.createFolder(name);
			}
			folders.put(id, dataFolder);
			folderTypes.put(dataFolder, folderType);
		}
	}

	private void createDataBeans() {
		
		for (DataType dataType : this.sessionType.getData()) {
			String name = dataType.getName();
			String id = dataType.getId();
						
			// check for unique session id
			if (getDataItem(id) != null) {
				logger.warn("duplicate data bean id: " + id + " , ignoring data bean: " + name);
				continue;
			}
			
			// check that data id exists
			String dataId = dataType.getDataId();
			if (dataId == null) {
				logger.warn("could not load data bean: " + name + " due to missing data id");
				throw new RuntimeException("trying to load data without data id");
			} 
									
			DataBean dataBean;
			try {
				
				/* Don't ask content length from filebroker at this point,
				 * but do it later in parallel along the type tags.
				 */				
				dataBean = dataManager.createDataBean(name, dataId, false);		
				
				for (LocationType location : dataType.getLocation()) {
					
					String urlString = location.getUrl();
					URL url = null;
					try {
						url = new URL(urlString);
					} catch (MalformedURLException e) {
						logger.warn("could not parse url: "  + urlString + " for data bean: " + name);
						continue;
					}

					if (StorageMethod.LOCAL_SESSION_ZIP.toString().equals(location.getMethod()) && !isDatalessSession) {
						// data is inside the session file, use the url for the real session file 
						url = new URL(sessionFile.toURI().toURL(), "#" + url.getRef());
					}

					dataManager.addContentLocationForDataBean(dataBean, StorageMethod.valueOfConverted(location.getMethod()), url);
				}
				
				// Set file size from metadata. If there are external
				// ContentLocations, the size must match.
				dataManager.setOrVerifyContentLength(dataBean, dataType.getSize());
				// set checksum from the metadata, but the checksum of the real file is calculated only 
				// later during possible network transfers
				dataManager.setOrVerifyChecksum(dataBean, dataType.getChecksum());
				
				Integer x = dataType.getLayoutX();
				Integer y = dataType.getLayoutY();
				
				if (x != null && y != null) {
					if (xOffset != null) {
						x += xOffset;
					}
					dataBean.setPosition(x, y);
				}
			
			} catch (Exception e) {
				Session.getSession().getApplication().reportExceptionThreadSafely(new Exception("error while opening file " + name, e));
				logger.warn("could not create data bean: " + name);
				continue;
			}
			
			// creation time
			XMLGregorianCalendar xmlCalendar = dataType.getCreationTime();
			if (xmlCalendar != null) {
				dataBean.setCreationDate(xmlCalendar.toGregorianCalendar().getTime());
			}

			dataBean.setNotes(dataType.getNotes());
			dataBean.setContentType(dataManager.guessContentType(dataBean.getName()));
			
			dataBeans.put(id, dataBean);
			dataTypes.put(dataBean, dataType);
		}		
	}

	
	private void createOperations() {
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

			// module
			operationRecord.setModule(operationType.getModule());
			
			// parameters
			for (ParameterType parameterType : operationType.getParameter()) {
				operationRecord.addParameter(parameterType.getName().getId(), parameterType.getName().getDisplayName(), parameterType.getName().getDescription(), parameterType.getValue());
			}

			// source code
			String sourceCodeFileName = operationType.getSourceCodeFile();
			if (sourceCodeFileName != null && !sourceCodeFileName.isEmpty()) {
				String sourceCode = null;
				try {
					sourceCode = getSourceCode(sourceCodeFileName);
				} catch (Exception e) {
					logger.warn("could not load source code from " + sourceCodeFileName);
				}

				// might be null, is ok
				operationRecord.setSourceCode(sourceCode);
			}

			// update names, category from the current version of the tool
			ClientApplication application = Session.getSession().getApplication();
			OperationDefinition currentTool = null;
			
			if (application != null) { //there is no client application when filebroker handles example sessions				
				currentTool = application.getOperationDefinitionBestMatch(operationRecord.getNameID().getID(), operationRecord.getModule(), operationRecord.getCategoryName());
			}
				
			if (currentTool != null) {
				operationRecord.getNameID().setDisplayName(currentTool.getDisplayName());
				operationRecord.getNameID().setDescription(currentTool.getDescription());
				if (currentTool.getCategory().getModule() != null) {
					operationRecord.setModule(currentTool.getCategory().getModule().getModuleName());
				}
				operationRecord.setCategoryName(currentTool.getCategoryName());
				operationRecord.setCategoryColor(currentTool.getCategory().getColor());

				for (ParameterRecord parameterRecord : operationRecord.getParameters()) {
					Parameter currentParameter = currentTool.getParameter(parameterRecord.getNameID().getID());
					if (currentParameter != null) {
						parameterRecord.getNameID().setDisplayName(currentParameter.getDisplayName());
						parameterRecord.getNameID().setDescription(currentParameter.getDescription());
					}
				}
			}

			// job id for continuation
			operationRecord.setJobId(operationType.getJobId());

			// start and end times
			XMLGregorianCalendar startTimeXML = operationType.getStartTime();
			if (startTimeXML != null) {
				operationRecord.setStartTime(startTimeXML.toGregorianCalendar().getTime());
			}
			XMLGregorianCalendar endTimeXML = operationType.getStartTime();
			if (startTimeXML != null) {
				operationRecord.setEndTime(endTimeXML.toGregorianCalendar().getTime());
			}
			
			// store the operation record
			operationRecords.put(operationSessionId, operationRecord);
			operationTypes.put(operationRecord, operationType);
		}
	}

	
	private void linkDataItemChildren(DataFolder parent) {
		
		ArrayList<DataItem> children = new ArrayList<>();
		ArrayList<DataFolder> folders = new ArrayList<>();
		
		for (String childId : folderTypes.get(parent).getChild()) {
			
			// check that the referenced data item exists
			DataItem child = getDataItem(childId);
			if (child == null) {
				logger.warn("child with id: " + childId + " does not exist");
				continue;
			}

			// add as a child
			children.add(child);
			
			// recursively go inside folders
			if (child instanceof DataFolder) {
				folders.add((DataFolder) child);
			}
		}
		
		// connect children in parallel
		dataManager.connectChildren(children, parent);
		
		for (DataFolder folder : folders) {
			linkDataItemChildren(folder);
		}
	}
	
	/**
	 * Link OperationRecords and DataBeans by adding real input DataBean references
	 * to OperationRecords.
	 * 
	 */
	private void linkInputsToOperations() {
		
		for (OperationRecord operationRecord : operationRecords.values()) {
			
			// get data bean ids from session data
			for (InputType inputType : operationTypes.get(operationRecord).getInput()) {

				String inputID = inputType.getData();
				
				// data bean exists
				
				if (inputID != null) {
					DataBean inputBean = dataBeans.get(inputID);
					
					// skip phenodata, it is bound automatically
					if (inputBean.queryFeatures("/phenodata/").exists()) {
						continue; 
					}
					// add the reference to the operation record
					operationRecord.addInput(createNameID(inputType.getName()), inputBean);
				}
				
				// data bean does not exist
				else {

					// try to skip phenodata, not reliable
					if (inputType.getName().getId().equals("phenodata.tsv")) {
						continue;
					}
					
					// add the reference to the operation record
					String dataId = inputType.getDataId();
					if (dataId != null) {
						operationRecord.addInput(createNameID(inputType.getName()), dataId);
					}
				}
			}
		}
	}

	
	/**
	 * Add links form DataBeans to the OperationRecord which created the DataBean.
	 * 
	 * If OperationRecord is not found, use unknown OperationRecord.
	 * 
	 */
	private void linkOperationsToOutputs() {
		
		for (DataBean dataBean : dataBeans.values()) {
			String operationId = dataTypes.get(dataBean).getResultOf();
			OperationRecord operationRecord = null;
			if (operationId != null) {
				operationRecord = operationRecords.get(operationId);
			}

			// if operation record is not found use dummy
			if (operationRecord == null) {
				operationRecord = OperationRecord.getUnkownOperationRecord();
			}
			
			dataBean.setOperationRecord(operationRecord);
		}
	}

	/**
	 * 
	 */
	private void linkDataBeans() {
		for (DataBean dataBean : dataBeans.values()) {
			for (LinkType linkType : dataTypes.get(dataBean).getLink()) {
				// if something goes wrong for this link, continue with others
				try {
					String targetID = linkType.getTarget();
					if (targetID != null) {
						DataBean targetBean = dataBeans.get(targetID);
						if (targetBean != null) {
							dataBean.addLink(Link.valueOf(linkType.getType()), targetBean);
						}
					}
					
				} catch (Exception e) {
					logger.warn("could not add link", e);
					continue;
				}
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

	private NameID createNameID(NameType name) {
		return new NameID(name.getId(), name.getDisplayName(), name.getDescription());
	}

	private String getSourceCode(String sourceCodeFileName) throws Exception {
		InputStream sourceCodeInputStream = null;
		StringWriter stringWriter = null;
		try {			
			sourceCodeInputStream = getStreamOfZipEntry(sourceCodeFileName);

			stringWriter = new StringWriter();
			IOUtils.copy(sourceCodeInputStream, new WriterOutputStream(stringWriter));
			stringWriter.flush();
		}
		finally {
			IOUtils.closeIfPossible(zipFile);
			IOUtils.closeIfPossible(zipStream);
			IOUtils.closeIfPossible(stringWriter);
		}

		return stringWriter.toString();
	}

	public List<OperationRecord> loadSession() throws Exception {
		
		// parse metadata to jaxb classes
		parseMetadata();

		// create the basic objects from the jaxb classes 
		createFolders();
		createDataBeans();
		createOperations();
		linkOperationsToOutputs();
				
		linkDataItemChildren(dataManager.getRootFolder());
		linkDataBeans();
		linkInputsToOperations();
		
		this.sessionNotes = sessionType.getNotes();
		return getUnfinishedOperations();
	}
	
	public void setXOffset(Integer xOffset) {
		this.xOffset = xOffset;
	}
	
	public String getSessionNotes() {
		return this.sessionNotes;
	}

	public List<OperationRecord> getUnfinishedOperations() {
		
		ArrayList<OperationRecord> unfinished = new ArrayList<>();
		
		for (OperationRecord operationRecord : this.operationRecords.values()) {
			String jobId = operationRecord.getJobId();
			if (jobId != null) {
				unfinished.add(operationRecord);
			}
		}
		return unfinished;
	}
}
